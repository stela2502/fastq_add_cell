use clap::Parser;
use needletail::parse_fastx_reader;
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Cli {
    #[arg(long, short)]
    cell: PathBuf,

    #[arg(long, short = '1')]
    r1: PathBuf,

    #[arg(long, short = '2')]
    r2: Option<PathBuf>,

    #[arg(long)]
    from_char: Option<usize>,

    #[arg(long)]
    to_char: Option<usize>,

    #[arg(long)]
    recomp: bool,
}

fn check_pigz_installed() -> Result<(), String> {
    match Command::new("pigz").arg("--version").output() {
        Ok(output) if output.status.success() => Ok(()),
        _ => Err("Error: `pigz` is not installed or not found in PATH.".to_string()),
    }
}


fn process_cell_sequence(
    seq: &[u8],
    from_char: Option<usize>,
    to_char: Option<usize>,
    revcomp: bool,
) -> String {
    let start = from_char.unwrap_or(0);
    let end = to_char.unwrap_or(seq.len());

    let clipped = if start < end && end <= seq.len() {
        &seq[start..end]
    } else {
        &seq
    };

    let final_seq = if revcomp {
        clipped
            .iter()
            .rev()
            .map(|b| match b {
                b'A' | b'a' => b'T',
                b'T' | b't' => b'A',
                b'G' | b'g' => b'C',
                b'C' | b'c' => b'G',
                b'N' | b'n' => b'N',
                _ => b'N',
            })
            .collect::<Vec<u8>>()
    } else {
        clipped.to_vec()
    };

    String::from_utf8_lossy(&final_seq).to_string()
}

fn fastq_record_to_string(id: &[u8], desc: &[u8], seq: &[u8], qual: Option<&[u8]>) -> String {
    let mut s = String::new();
    s.push('@');
    let cleaned: String = id.iter().map(|&b| if b == b' ' { b':' } else { b })
        .collect::<Vec<u8>>()
        .into_iter()
        .map(|b| b as char)
        .collect();
    s.push_str(&cleaned);
    s.push(':');
    s.push_str(&String::from_utf8_lossy(desc));
    s.push('\n');
    s.push_str(&String::from_utf8_lossy(seq));
    s.push('\n');
    s.push('+');
    s.push('\n');
    if let Some(q) = qual {
        s.push_str(&String::from_utf8_lossy(q));
    }
    s.push('\n');
    s
}

fn start_pigz_reader(path: &Path) -> std::io::Result<BufReader<std::process::ChildStdout>> {
    let child = Command::new("pigz")
        .arg("-dc")
        .arg(path)
        .stdout(Stdio::piped())
        .spawn()?;
    Ok(BufReader::new(child.stdout.expect("No stdout from pigz")))
}

fn start_pigz_writer(output_path: &Path, threads: usize) -> std::io::Result<(std::process::Child, BufWriter<std::process::ChildStdin>)> {
    let out_file = File::create(output_path)?;
    let mut child = Command::new("pigz")
        .arg("-p").arg(threads.to_string())
        .arg("-c")
        .stdin(Stdio::piped())
        .stdout(Stdio::from(out_file))
        .spawn()?;
    let child_stdin = child.stdin
        .take()
        .ok_or_else(|| std::io::Error::new(std::io::ErrorKind::Other, "Failed to capture pigz stdin"))?;

    Ok((child, BufWriter::new(child_stdin)))
}

fn make_output_name(input_path: &Path) -> PathBuf {
    let filename = input_path.file_name().unwrap().to_string_lossy();
    if filename.ends_with(".fastq.gz") {
        let base = filename.trim_end_matches(".fastq.gz");
        input_path.with_file_name(format!("{}_cells_added.fastq.gz", base))
    } else {
        input_path.with_file_name(format!("{}_cells_added.fastq.gz", filename))
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli = Cli::parse();

    check_pigz_installed()?;

    let cell_reader = start_pigz_reader(&cli.cell)?;
    let r1_reader = start_pigz_reader(&cli.r1)?;
    let r2_reader = if let Some(ref r2) = cli.r2 {
        Some(start_pigz_reader(r2)?)
    } else {
        None
    };

    let mut cell_parser = parse_fastx_reader(cell_reader)?;
    let mut r1_parser = parse_fastx_reader(r1_reader)?;
    let mut r2_parser = if let Some(r) = r2_reader {
        Some(parse_fastx_reader(r)?)
    } else {
        None
    };

    let (mut r1_child, mut r1_writer) = start_pigz_writer(&make_output_name(&cli.r1), 4)?;
    let (r2_child, mut r2_writer) = if let Some(ref r2) = cli.r2 {
        let (child, writer) = start_pigz_writer(&make_output_name(r2), 4)?;
        (Some(child) ,Some(writer))
    } else {
        (None, None)
    };

    loop {
        let cell_rec = match cell_parser.next() {
            Some(Ok(r)) => r,
            _ => break,
        };
        let r1_rec = match r1_parser.next() {
            Some(Ok(r)) => r,
            _ => break,
        };
        let r2_rec = match r2_parser.as_mut().and_then(|p| p.next()) {
            Some(Ok(r)) => Some(r),
            _ => None,
        };

        let cell_seq = process_cell_sequence(&cell_rec.seq(), cli.from_char, cli.to_char, cli.recomp);
        let cell_bytes = cell_seq.as_bytes();

        let r1_str = fastq_record_to_string(&r1_rec.id(), cell_bytes, &r1_rec.seq(), r1_rec.qual());
        r1_writer.write_all(r1_str.as_bytes())?;

        if let (Some(writer), Some(r2)) = (&mut r2_writer, r2_rec) {
            let r2_str = fastq_record_to_string(&r2.id(), cell_bytes, &r2.seq(), r2.qual());
            writer.write_all(r2_str.as_bytes())?;
        }
    }

    drop(r1_writer);
    if let Some(writer) = r2_writer { drop(writer); }

    r1_child.wait()?;
    if let Some(mut child) = r2_child { child.wait()?; }

    Ok(())
}
