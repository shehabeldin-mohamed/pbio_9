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


def format_fasta(seq_id: str, description: str, sequence: str, line_width: int = 80) -> str:
    """Returns a formatted FASTA record as a string."""
    # Format the header
    header = f">{seq_id}"
    if description:
        header += f" {description}"

    lines = [header]

    # Break the sequence into lines of fixed width
    for i in range(0, len(sequence), line_width):
        lines.append(sequence[i:i + line_width])

    lines.append("")
    lines.append("# EOF_1")

    return "\n".join(lines)


def main():
    """Main function handling user interaction, generation, and file output."""
    # 1. Gather sequence length
    length = validate_positive_int("Enter sequence length: ")

    # 2. Gather and validate Sequence ID
    while True:
        seq_id = input("Enter sequence ID: ")
        if " " in seq_id or "\t" in seq_id:
            print("Error: ID cannot contain whitespace. Please try again.")
        else:
            break

    # 3. Gather optional description and name
    description = input("Enter a description of the sequence: ")
    name = input("Enter your name: ")

    # 4. Generate the biological sequence
    pure_sequence = generate_sequence(length)

    # 5. Calculate statistics BEFORE inserting the name
    stats = calculate_stats(pure_sequence)

    # 6. Insert the name visually
    final_sequence = insert_name(pure_sequence, name)

    # 7. Format and save to file
    fasta_output = format_fasta(seq_id, description, final_sequence)
    filename = f"{seq_id}.fasta"

    try:
        with open(filename, "w", encoding="utf-8") as file:
            file.write(fasta_output)
        print(f"\nSequence saved to file: {filename}\n")
    except IOError as e:
        print(f"\nError saving file: {e}")
        sys.exit(1)

    # 8. Print statistics and sample output
    print(f"Sequence statistics (n={length}):")
    print(f"A: {stats['A']:.2f}%")
    print(f"C: {stats['C']:.2f}%")
    print(f"G: {stats['G']:.2f}%")
    print(f"T: {stats['T']:.2f}%")
    print(f"GC-content: {stats['gc_ratio_A']:.2f}%\n")

    print(f"Sample contents of the {filename} file:")
    print(fasta_output)


if __name__ == "__main__":
    main()