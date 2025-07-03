use clap::Parser;
use needletail::{parse_fastx_file};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use flate2::write::GzEncoder;
use flate2::Compression;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Cli {
    /// FASTQ file for cell (required)
    #[arg(long, short)]
    cell: PathBuf,

    /// FASTQ file for R1 (required)
    #[arg(long, short = '1')]
    r1: PathBuf,

    /// FASTQ file for R2 (optional)
    #[arg(long, short = '2')]
    r2: Option<PathBuf>,

    /// Ommit some nulceotides from the start of the cell read?
    #[arg(long)]
    from_char: Option<usize>,

    /// Stop reading from the cell read prematurely?
    #[arg(long)]
    to_char: Option<usize>,

    /// Does the Cell read need to be reverse complemented?
    #[arg(long)]
    recomp: bool
}

pub fn process_cell_sequence(
    seq: &[u8],
    from_char: Option<usize>,
    to_char: Option<usize>,
    revcomp: bool,
) -> String {
    // Apply front and back clipping
    let start = from_char.unwrap_or(0);
    let end = to_char.unwrap_or(seq.len());

    // Bounds check
    let clipped = if start < end && end <= seq.len() {
        &seq[start..end]
    } else {
        &seq
    };

    // Reverse complement if needed
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
                _ => b'N', // fallback for unexpected characters
            })
            .collect::<Vec<u8>>()
    } else {
        clipped.to_vec()
    };

    // Convert to String
    String::from_utf8_lossy(&final_seq).to_string()
}

fn fastq_record_to_string(
    id: &[u8],
    desc: &[u8],
    seq: &[u8],
    qual: Option<&[u8]>
) -> String {
    let mut s = String::new();

    s.push('@');
    s.push_str(&String::from_utf8_lossy(id));
    s.push(' ');
    s.push_str(&String::from_utf8_lossy(desc));
    s.push('\n');

    s.push_str(&String::from_utf8_lossy(seq));
    s.push('\n');

    s.push('+');
    s.push('\n');
    if let Some(q) =  qual {
        s.push_str(&String::from_utf8_lossy(q));
    }
    s.push('\n');

    s
}

fn make_output_name(input_path: &Path) -> PathBuf {
    let filename = input_path.file_name().unwrap().to_string_lossy();

    // Handle .gz compressed files
    if filename.ends_with(".fastq.gz") {
        let base = filename.trim_end_matches(".fastq.gz");
        input_path.with_file_name(format!("{}_cells_added.fastq.gz", base))
    } else if filename.ends_with(".fq.gz") {
        let base = filename.trim_end_matches(".fq.gz");
        input_path.with_file_name(format!("{}_cells_added.fq.gz", base))
    }
    // Handle uncompressed files
    else if filename.ends_with(".fastq") {
        let base = filename.trim_end_matches(".fastq");
        input_path.with_file_name(format!("{}_cells_added.fastq", base))
    } else if filename.ends_with(".fq") {
        let base = filename.trim_end_matches(".fq");
        input_path.with_file_name(format!("{}_cells_added.fq", base))
    }
    // Fallback: just append _cells_added if no recognized extension
    else {
        let stem = input_path.file_stem().unwrap().to_string_lossy();
        let ext = input_path.extension().map(|e| e.to_string_lossy()).unwrap_or_else(|| "".into());
        if ext.is_empty() {
            input_path.with_file_name(format!("{}_cells_added", stem))
        } else {
            input_path.with_file_name(format!("{}_cells_added.{}", stem, ext))
        }
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli = Cli::parse();

    // Open cell FASTQ reader
    let mut cell_reader = parse_fastx_file(&cli.cell).expect("failed to read from cell fastq");
    let mut r1_reader = parse_fastx_file(&cli.r1).expect("failed to read from r1 fastq");

    // Optional r2 reader
    let mut r2_reader = match &cli.r2 {
        Some(p) => Some( parse_fastx_file(p).expect("failed to read from r2 fastq") ),
        None => None,
    };

    // Create output writers with gzip compression
    //let out_cell = File::create(make_output_name(&cli.cell))?;
    let out_r1 = File::create(make_output_name(&cli.r1))?;
    let out_r2 = cli.r2.as_ref().map(|p| File::create(make_output_name(p))).transpose()?;

    //let mut cell_writer = FastqWriter::new(GzEncoder::new(BufWriter::new(out_cell), Compression::default()));
    let buf_r1_writer = BufWriter::new(out_r1);
    let mut r1_writer = GzEncoder::new(buf_r1_writer, Compression::default());

    let mut r2_writer: Option<GzEncoder<BufWriter<File>>> = match out_r2 {
        Some(file) => {
            let buf_r2_writer = BufWriter::new(file);
            let r2_w = GzEncoder::new(buf_r2_writer, Compression::default());
            Some(r2_w)
        },
        None => None,
    };

    loop {
        // Read one record from each input; break if any run out (except r2 which is optional)
        let cell_rec = match cell_reader.next() {
            Some(Ok(r)) => r,
            _ => break,
        };
        let r1_rec = match r1_reader.next() {
            Some(Ok(r)) => r,
            _ => break,
        };
        let r2_rec = match &mut r2_reader {
            Some(r2r) => match r2r.next() {
                Some(Ok(r)) => Some(r),
                _ => break,
            },
            None => None,
        };

        // Get the cell sequence string
        let cell_seq = process_cell_sequence(&cell_rec.seq(), cli.from_char, cli.to_char, cli.recomp );
        let cell = cell_seq.as_bytes();

        // Write cell read unchanged to output
        //cell_writer.write(cell_rec.id(), cell_rec.desc(), cell_rec.seq(), cell_rec.qual())?;

        // Modify r1 ID by appending |cell:<cell_seq>
        let r1_str = fastq_record_to_string( &r1_rec.id(), &cell, &r1_rec.seq(), r1_rec.qual() );
        r1_writer.write( r1_str.as_bytes() )?;

        // Modify and write r2 if present
        if let (Some(r2_writer), Some(r2_rec)) = (r2_writer.as_mut(), r2_rec) {
            let r2_str = fastq_record_to_string( &r2_rec.id(), &cell, &r2_rec.seq(), r2_rec.qual() );
            r2_writer.write(r2_str.as_bytes())?;
        }
    }

    Ok(())
}



#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_process_cell_sequence_clipping_and_revcomp() {
        let seq = b"ACGTTAGC";

        // Clip 2 from front and 2 from back => "GTT"
        // Revcomp => "AAC"
        let result = process_cell_sequence(seq, Some(2), Some(5), true);
        assert_eq!(result, "AAC", "Expected reverse complement of 'GTT' to be 'CTT'");
    }

    #[test]
    fn test_no_clip_no_revcomp() {
        let seq = b"ATGC";
        let result = process_cell_sequence(seq, None, None, false);
        assert_eq!(result, "ATGC", "Should return original sequence unchanged");
    }

    #[test]
    fn test_only_revcomp() {
        let seq = b"ATGC";
        let result = process_cell_sequence(seq, None, None, true);
        assert_eq!(result, "GCAT", "Expected reverse complement of 'ATGC' to be 'GCAT'");
    }

    #[test]
    fn test_clip_out_of_bounds() {
        let seq = b"ATGC";
        let result = process_cell_sequence(seq, Some(10), Some(0), false);
        assert_eq!(result, "ATGC", "Clipping beyond bounds should return original");
    }

    #[test]
    fn test_clip_to_single_base_and_revcomp() {
        let seq = b"CCACCC";
        let result = process_cell_sequence(seq, Some(2), Some(3), true); // should be "A", revcomp -> "T"
        assert_eq!(result, "T", "Expected reverse complement of 'A' to be 'T'");
    }
}