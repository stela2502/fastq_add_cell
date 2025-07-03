use std::process::Command;
use std::fs::File;
use std::io::{BufReader};
use flate2::read::MultiGzDecoder;
use needletail::parse_fastx_reader;

#[test]
fn test_fastq_add_cell_ids() {
    let cell_path = "testData/Cell.fastq.gz";
    let r1_path = "testData/R1.fastq.gz";
    let r2_path = "testData/R2.fastq.gz";

    // Run the binary
    let output = Command::new("target/release/fastq_add_cell")
        .args(&["-c", cell_path, "-1", r1_path, "-2", r2_path])
        .output()
        .expect("failed to execute fastq_add_cell");

    assert!(output.status.success());

    // Determine output file names (assuming naming logic)
    let out_r1 = "testData/R1_cells_added.fastq.gz";
    let out_r2 = "testData/R2_cells_added.fastq.gz";

    // Read all cell sequences into a vector
    let mut cell_sequences = Vec::new();
    {
        let file = File::open(cell_path).unwrap();
        let reader = BufReader::new(MultiGzDecoder::new(file));
        let mut fastx = parse_fastx_reader(reader).unwrap();

        while let Some(record) = fastx.next() {
            let seqrec = record.expect("valid cell FASTQ record");
            let seq = std::str::from_utf8(&seqrec.seq()).unwrap().to_string();
            cell_sequences.push(seq);
        }
    }

    // Check that each R1 read ID ends with the corresponding cell sequence
    {

        let file = File::open(out_r1).unwrap();
        let reader = BufReader::new(MultiGzDecoder::new(file));
        let mut fastx = parse_fastx_reader(reader).unwrap();

        let mut i = 0;
        while let Some(record) = fastx.next() {
            let rec = record.expect("valid R1 FASTQ record");
            let id = std::str::from_utf8(rec.id()).unwrap();
            let cell = &cell_sequences[i];
            assert!(
                id.contains(cell),
                "R1 ID does not contain expected cell barcode: id='{}', cell='{}'",
                id,
                cell
            );
            i += 1;
        }
    }

    // Check that each R2 read ID ends with the corresponding cell sequence
    {
        let file = File::open(out_r2).unwrap();
        let reader = BufReader::new(MultiGzDecoder::new(file));
        let mut fastx = parse_fastx_reader(reader).unwrap();

        let mut i = 0;
        while let Some(record) = fastx.next() {
            let rec = record.expect("valid R2 FASTQ record");
            let id = std::str::from_utf8(rec.id()).unwrap();
            let cell = &cell_sequences[i];
            assert!(
                id.contains(cell),
                "R2 ID does not contain expected cell barcode: id='{}', cell='{}'",
                id,
                cell
            );
            i +=1;
        }
    }
}
