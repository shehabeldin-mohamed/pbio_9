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
    # Use only uppercase letters to ensure inserted names don't skew the math
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

    # specific validator for the end of the file
    if include_eof:
        lines.append("")
        lines.append("# EOF_1")

    return "\n".join(lines)


def generate_batch(num_sequences: int, length: int, base_id: str, description: str, name: str) -> tuple[str, list]:
    """
    Generates multiple sequences in a loop, appends them into a single
    multi-FASTA string, and returns the string alongside a list of statistics.
    """
    multi_fasta_blocks = []
    batch_stats = []

    for i in range(1, num_sequences + 1):
        # Create unique ID padded with zeros (e.g., Seq_001, Seq_002)
        current_id = f"{base_id}_{i:03d}"

        pure_seq = generate_sequence(length)
        stats = calculate_stats(pure_seq)
        batch_stats.append((current_id, stats))

        final_seq = insert_name(pure_seq, name)

        # Only attach the required EOF marker to the final sequence in the batch loop
        is_last = (i == num_sequences)
        fasta_block = format_fasta(current_id, description, final_seq, include_eof=is_last)

        multi_fasta_blocks.append(fasta_block)

    # Join distinct FASTA blocks with a newline separator
    return "\n".join(multi_fasta_blocks), batch_stats


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

    # Process sequence generation based on selected mode
    filename = f"{base_id}.fasta" if mode == 1 else f"{base_id}_batch.fasta"

    if mode == 1:
        pure_sequence = generate_sequence(length)
        stats = calculate_stats(pure_sequence)
        final_sequence = insert_name(pure_sequence, name)
        fasta_output = format_fasta(base_id, description, final_sequence)
        all_stats = [(base_id, stats)]
    else:
        fasta_output, all_stats = generate_batch(num_sequences, length, base_id, description, name)

    try:
        with open(filename, "w", encoding="utf-8") as file:
            file.write(fasta_output)
        print(f"\nSequence(s) saved to file: {filename}\n")
    except IOError as e:
        print(f"\nError saving file: {e}")
        sys.exit(1)

    # Print statistics for first sequence only
    print(f"Sequence statistics calculated for {len(all_stats)} sequence(s).")
    print(f"Showing stats for the first sequence ({all_stats[0][0]}, n={length}):")
    first_stats = all_stats[0][1]
    print(f"A: {first_stats['A']:.2f}%")
    print(f"C: {first_stats['C']:.2f}%")
    print(f"G: {first_stats['G']:.2f}%")
    print(f"T: {first_stats['T']:.2f}%")
    print(f"GC-content: {first_stats['gc_ratio_A']:.2f}%\n")

    # Print sample file contents if output is less than 600 otherwise truncated
    print(f"Sample contents of the {filename} file:")
    if len(fasta_output) > 600:
        print(fasta_output[:600] + "\n... [TRUNCATED] ...")
    else:
        print(fasta_output)


if __name__ == "__main__":
    main()