import random


def generate_dna_sequence(length):
    nucleotides = ['A', 'C', 'G', 'T']
    return ''.join(random.choice(nucleotides) for _ in range(length))


def insert_name_in_sequence(sequence, name):
    position = random.randint(0, len(sequence))
    return sequence[:position] + name + sequence[position:]


def calculate_statistics(sequence):
    length = len(sequence)
    a_count = sequence.count('A')
    c_count = sequence.count('C')
    g_count = sequence.count('G')
    t_count = sequence.count('T')

    a_percentage = (a_count / length) * 100
    c_percentage = (c_count / length) * 100
    g_percentage = (g_count / length) * 100
    t_percentage = (t_count / length) * 100

    c_g_ratio = (c_count + g_count) / (a_count + t_count) if (a_count + t_count ) != 0 else 0

    stats = {
        'A': a_percentage,
        'C': c_percentage,
        'G': g_percentage,
        'T': t_percentage,
        'C/G to A/T ratio': c_g_ratio
    }
    return stats


def save_to_fasta(sequence_id, description, sequence):
    file_name = f"{sequence_id}.fasta"
    with open(file_name, "w") as fasta_file:
        fasta_file.write(f">{sequence_id} {description}\n")
        fasta_file.write(sequence + "\n")
    print(f"Plik {file_name} zapisany.")


def main():
    #ORIGINAL:
    #length = int(input("Podaj długość sekwencji: "))

    #MODIFIED ()
    MAX_SEQUENCE_LENGTH = 10**6
    while True:
        try:
            length = int(input("Podaj długość sekwencji: "))
            if length <= 0:
                raise ValueError("Długość musi być liczbą dodatnią.")
            if length > MAX_SEQUENCE_LENGTH:
                raise ValueError(f"Długość sekwencji nie może przekroczyć {MAX_SEQUENCE_LENGTH}.")
            break
        except ValueError as e:
            print(f"Błąd: {e}. Spróbuj ponownie.")

    sequence_id = input("Podaj nazwę (ID) sekwencji: ")
    description = input("Podaj opis sekwencji: ")
    name = input("Podaj imię: ")

    dna_sequence = generate_dna_sequence(length)

    sequence_with_name = insert_name_in_sequence(dna_sequence, name)

    stats = calculate_statistics(dna_sequence)

    # ORIGINAL:
    # print("\nStatystyki sekwencji DNA:")
    # for key, value in stats.items():
    #     if key == 'C/G to A/T ratio':
    #         print(f"%CG: {value:.2f}")
    #     else:
    #         print(f"{key}: {value:.2f}%")

    #MODIFIED
    print("\nStatystyki sekwencji DNA:")
    file_name_statistics = f"{sequence_id}_statistics.txt"
    try:
        with open(file_name_statistics, "w") as stats_file:
            for key, value in stats.items():
                if key == 'C/G to A/T ratio':
                    stats_file.write(f"%CG: {value:.2f}\n")
                    print(f"%CG: {value:.2f}")
                else:
                    stats_file.write(f"{key}: {value:.2f}%\n")
                    print(f"{key}: {value:.2f}%")
        print(f"Statystyki zapisane do pliku {file_name_statistics}.")
    except IOError as e:
        print(f"Nie udało się zapisać statystyk: {e}")


    save_to_fasta(sequence_id, description, sequence_with_name)



if __name__ == "__main__":
    main()
