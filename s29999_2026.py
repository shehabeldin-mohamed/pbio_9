# Album number: s29999
# Date: 01-05-2026
# Description: Generates random DNA sequences, calculates stats, finds motifs,
# generates complements/translations, and saves everything to a FASTA file.

import random
import sys

def validate_positive_int(prompt: str, min_val: int = 1, max_val: int = 100_000) -> int:
    """
    Gets an integer from the user within a specified range.
    In case of an error, repeats the question.
    """
    while True:
        try:
            user_input = input(prompt)
            value = int(user_input)
            if min_val <= value <= max_val:
                return value
            else:
                print(f"Error: value must be an integer in the range [{min_val}, {max_val}].")
        except ValueError:
            print(f"Error: value must be an integer in the range [{min_val}, {max_val}].")


def generate_sequence(length: int) -> str:
    """Returns a random DNA sequence of the specified length."""
    nucleotides = ['A', 'C', 'G', 'T']
    return ''.join(random.choices(nucleotides, k=length))


def calculate_stats(sequence: str) -> dict:
    """
    Returns a dictionary of sequence statistics.
    Keys: "A", "C", "G", "T" (float values, %),
          "gc_ratio_A" (float value, %).
    """
    pure_seq = [base for base in sequence if base.isupper()]
    seq_len = len(pure_seq)

    if seq_len == 0:
        return {"A": 0.0, "C": 0.0, "G": 0.0, "T": 0.0, "gc_ratio_A": 0.0}

    stats = {
        "A": (pure_seq.count('A') / seq_len) * 100,
        "C": (pure_seq.count('C') / seq_len) * 100,
        "G": (pure_seq.count('G') / seq_len) * 100,
        "T": (pure_seq.count('T') / seq_len) * 100,
    }

    stats["gc_ratio_A"] = stats["G"] + stats["C"]

    return stats


def insert_name(sequence: str, name: str) -> str:
    """Inserts a name at a random position in the sequence in lowercase letters."""
    if not sequence:
        return name.lower()

    insert_pos = random.randint(0, len(sequence))
    return sequence[:insert_pos] + name.lower() + sequence[insert_pos:]


def format_fasta(seq_id: str, description: str, sequence: str, line_width: int = 80, include_eof: bool = True) -> str:
    """Returns a formatted FASTA record as a string."""
    header = f">{seq_id}"
    if description:
        header += f" {description}"

    lines = [header]

    # Break the sequence into lines of fixed width
    for i in range(0, len(sequence), line_width):
        lines.append(sequence[i:i + line_width])

    if include_eof:
        lines.append("")
        lines.append("# EOF_1")

    return "\n".join(lines)


def find_motif(sequence: str, motif: str) -> list[int]:
    """
    Searches for a user-specified motif in the sequence.
    Returns a list of 1-based indexing positions of occurrences.
    """
    positions = []
    if not motif:
        return positions

    motif_len = len(motif)
    for i in range(len(sequence) - motif_len + 1):
        if sequence[i:i + motif_len] == motif:
            positions.append(i + 1)
    return positions


def get_complement(sequence: str) -> str:
    """Generates the complementary DNA strand (A->T, C->G)."""
    complement_map = str.maketrans('ACGT', 'TGCA')
    return sequence.translate(complement_map)


def get_reverse_complement(sequence: str) -> str:
    """Generates the reverse complementary DNA strand."""
    return get_complement(sequence)[::-1]


def translate_to_protein(sequence: str) -> str:
    """
    Translates a DNA sequence into an amino acid sequence using a standard codon table.
    Stop codons are represented by '*'.
    """
    codon_table = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
        'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
    }

    protein = []
    # Process the sequence in chunks of 3 (codons)
    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i + 3]
        protein.append(codon_table.get(codon, '?'))  # '?' for incomplete/unknown codons

    return "".join(protein)


def process_single_sequence(seq_id: str, description: str, pure_seq: str, name: str) -> list[tuple]:
    """
    Helper function to generate all records for a given base sequence.
    Returns a list of tuples containing (ID, description, sequence).
    """
    records = []

    # Base Sequence (with name injected)
    final_base_seq = insert_name(pure_seq, name)
    records.append((seq_id, description, final_base_seq))

    # Complement Sequence (pure)
    comp_seq = get_complement(pure_seq)
    records.append((f"{seq_id}_comp", f"Complement of {seq_id}", comp_seq))

    # Reverse Complement Sequence (pure)
    rev_comp_seq = get_reverse_complement(pure_seq)
    records.append((f"{seq_id}_revcomp", f"Reverse complement of {seq_id}", rev_comp_seq))

    # Protein Translation (pure)
    protein_seq = translate_to_protein(pure_seq)
    records.append((f"{seq_id}_prot", f"Translation of {seq_id}", protein_seq))

    return records


def generate_batch(num_sequences: int, length: int, base_id: str, description: str, name: str, motif: str) -> tuple[
    str, list, list]:
    """
    Handles generation of all biological records, formats them, and returns
    the multi-FASTA string, sequence statistics, and motif occurrences.
    """
    multi_fasta_blocks = []
    batch_stats = []
    batch_motifs = []

    all_records_to_format = []

    for i in range(1, num_sequences + 1):
        current_id = f"{base_id}_{i:03d}" if num_sequences > 1 else base_id

        # Core generation and analysis
        pure_seq = generate_sequence(length)
        stats = calculate_stats(pure_seq)
        motif_positions = find_motif(pure_seq, motif)

        batch_stats.append((current_id, stats))
        if motif:
            batch_motifs.append((current_id, motif_positions))

        # Get base, complement, rev_comp, and translation records
        seq_records = process_single_sequence(current_id, description, pure_seq, name)
        all_records_to_format.extend(seq_records)

    # Format all collected records into FASTA
    for index, record in enumerate(all_records_to_format):
        rec_id, rec_desc, rec_seq = record
        # Ensure ONLY the absolute last record triggers the `# EOF_1` flag
        is_last = (index == len(all_records_to_format) - 1)
        fasta_block = format_fasta(rec_id, rec_desc, rec_seq, include_eof=is_last)
        multi_fasta_blocks.append(fasta_block)

    return "\n".join(multi_fasta_blocks), batch_stats, batch_motifs


def main():
    """Main function handling user interaction, generation modes, and file output."""
    print("--- Sequence Generator ---")
    print("1. Single sequence mode")
    print("2. Batch mode (multi-FASTA)")
    mode = validate_positive_int("Select mode (1 or 2): ", 1, 2)

    length = validate_positive_int("Enter sequence length: ")

    num_sequences = 1
    if mode == 2:
        num_sequences = validate_positive_int("Enter number of sequences to generate: ", 1, 10_000)

    while True:
        base_id = input("Enter sequence base ID: ")
        if " " in base_id or "\t" in base_id:
            print("Error: ID cannot contain whitespace. Please try again.")
        else:
            break

    description = input("Enter a description of the sequence: ")
    name = input("Enter your name: ")
    motif = input("Enter a motif to search for (e.g., ATG) or press Enter to skip: ").strip().upper()

    # Process sequence generation
    filename = f"{base_id}.fasta" if mode == 1 else f"{base_id}_batch.fasta"

    fasta_output, all_stats, all_motifs = generate_batch(num_sequences, length, base_id, description, name, motif)

    try:
        with open(filename, "w", encoding="utf-8") as file:
            file.write(fasta_output)
        print(f"\n[SUCCESS] Sequence(s) and their variants saved to file: {filename}\n")
    except IOError as e:
        print(f"\n[ERROR] Could not save file: {e}")
        sys.exit(1)

    # Print statistics
    print(f"--- Statistics ({len(all_stats)} sequence(s) generated) ---")
    first_stats = all_stats[0][1]
    print(f"Stats for {all_stats[0][0]} (n={length}):")
    print(
        f"  A: {first_stats['A']:.2f}% | C: {first_stats['C']:.2f}% | G: {first_stats['G']:.2f}% | T: {first_stats['T']:.2f}%")
    print(f"  GC-content: {first_stats['gc_ratio_A']:.2f}%\n")

    # Print motif search results
    if motif:
        print(f"--- Motif Search Results for '{motif}' ---")
        for seq_id, positions in all_motifs:
            if positions:
                print(f"  {seq_id}: Found at positions {positions}")
            else:
                print(f"  {seq_id}: No matches found.")
        print()

    # Print sample file contents if output is less than 600 otherwise truncated
    print(f"--- Sample contents of {filename} ---")
    if len(fasta_output) > 600:
        print(fasta_output[:600] + "\n... [TRUNCATED] ...\n# EOF_1")
    else:
        print(fasta_output)


if __name__ == "__main__":
    main()