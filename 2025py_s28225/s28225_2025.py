import random  # Importuje moduł do generowania liczb losowych (służy do tworzenia losowych sekwencji DNA i wybierania losowych pozycji)
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt  # Importuje bibliotekę do tworzenia wykresów (matplotlib) (używana do wizualizacji statystyk nukleotydów)


def generate_dna_sequence(length): # Definiuje funkcję o nazwie 'generate_dna_sequence' przyjmującą jeden argument: 'length' (długość sekwencji)
    # Definiuje możliwe nukleotydy
    nucleotides = ['A', 'C', 'G', 'T'] # Tworzy listę zawierającą cztery podstawowe nukleotydy DNA
    # Zwraca losową sekwencję DNA o zadanej długości
    return ''.join(random.choice(nucleotides) for _ in range(length)) # Generuje ciąg znaków (sekwencję DNA) poprzez losowe wybieranie nukleotydów z listy 'nucleotides' 'length' razy i łączenie ich w jeden string


def insert_name_in_sequence(sequence, name): # Definiuje funkcję o nazwie 'insert_name_in_sequence' przyjmującą dwa argumenty: 'sequence' (sekwencja DNA) i 'name' (imię do wstawienia)
    # Wybiera losową pozycję w sekwencji
    position = random.randint(0, len(sequence)) # Generuje losową liczbę całkowitą od 0 do długości sekwencji (włącznie), która będzie pozycją wstawienia imienia
    # Wstawia imię w wybrane miejsce w sekwencji
    return sequence[:position] + name + sequence[position:] # Dzieli sekwencję na dwie części w miejscu 'position', wstawia 'name' pomiędzy nie i zwraca nową, połączoną sekwencję


def calculate_statistics(sequence): # Definiuje funkcję o nazwie 'calculate_statistics' przyjmującą jeden argument: 'sequence' (sekwencja DNA)
    length = len(sequence)  # Oblicza długość sekwencji i przypisuje ją do zmiennej 'length'
    a_count = sequence.count('A')  # Liczy wystąpienia nukleotydu 'A' w sekwencji i przypisuje wynik do 'a_count'
    c_count = sequence.count('C')  # Liczy wystąpienia nukleotydu 'C' w sekwencji i przypisuje wynik do 'c_count'
    g_count = sequence.count('G')  # Liczy wystąpienia nukleotydu 'G' w sekwencji i przypisuje wynik do 'g_count'
    t_count = sequence.count('T')  # Liczy wystąpienia nukleotydu 'T' w sekwencji i przypisuje wynik do 't_count'

    # Oblicza procentową zawartość każdego nukleotydu
    a_percentage = (a_count / length) * 100 if length > 0 else 0 # Oblicza procentową zawartość 'A', zabezpieczając przed dzieleniem przez zero, jeśli sekwencja jest pusta
    c_percentage = (c_count / length) * 100 if length > 0 else 0 # Oblicza procentową zawartość 'C', zabezpieczając przed dzieleniem przez zero
    g_percentage = (g_count / length) * 100 if length > 0 else 0 # Oblicza procentową zawartość 'G', zabezpieczając przed dzieleniem przez zero
    t_percentage = (t_count / length) * 100 if length > 0 else 0 # Oblicza procentową zawartość 'T', zabezpieczając przed dzieleniem przez zero

    # Oblicza stosunek zawartości C+G do A+T (jeśli A+T ≠ 0)
    c_g_ratio = (c_count + g_count) / (a_count + t_count) if (a_count + t_count) != 0 else 0 # Oblicza stosunek sumy C i G do sumy A i T, zabezpieczając przed dzieleniem przez zero, jeśli suma A i T wynosi 0

    # Zwraca słownik ze statystykami
    stats = { # Tworzy słownik o nazwie 'stats' przechowujący obliczone statystyki
        'A': a_percentage, # Klucz 'A' z wartością procentowej zawartości adeniny
        'C': c_percentage, # Klucz 'C' z wartością procentowej zawartości cytozyny
        'G': g_percentage, # Klucz 'G' z wartością procentowej zawartości guaniny
        'T': t_percentage, # Klucz 'T' z wartością procentowej zawartości tyminy
        'C/G to A/T ratio': c_g_ratio # Klucz 'C/G to A/T ratio' z wartością obliczonego stosunku GC do AT
    }
    return stats # Zwraca słownik 'stats'


def save_to_fasta(sequence_id, description, sequence): # Definiuje funkcję o nazwie 'save_to_fasta' przyjmującą trzy argumenty: 'sequence_id' (identyfikator sekwencji), 'description' (opis sekwencji) i 'sequence' (sekwencja DNA)
    # Tworzy nazwę pliku na podstawie ID
    file_name = f"{sequence_id}.fasta" # Tworzy nazwę pliku w formacie 'ID_sekwencji.fasta'
    # Otwiera plik FASTA do zapisu
    with open(file_name, "w") as fasta_file: # Otwiera plik o nazwie 'file_name' w trybie zapisu ("w") i przypisuje go do zmiennej 'fasta_file'. Użycie 'with' zapewnia automatyczne zamknięcie pliku.
        # Zapisuje nagłówek FASTA (ID i opis)
        fasta_file.write(f">{sequence_id} {description}\n") # Zapisuje linię nagłówka formatu FASTA, zaczynającą się od '>', po którym następuje ID i opis, zakończone znakiem nowej linii
        # Zapisuje sekwencję DNA z imieniem
        fasta_file.write(sequence + "\n") # Zapisuje właściwą sekwencję DNA do pliku, zakończoną znakiem nowej linii
    print(f"Plik {file_name} zapisany.")  # Informuje użytkownika o pomyślnym zapisaniu pliku, wyświetlając jego nazwę


def main(): # Definiuje główną funkcję programu o nazwie 'main'
    # Maksymalna dozwolona długość sekwencji
    MAX_SEQUENCE_LENGTH = 10**6 # Definiuje stałą 'MAX_SEQUENCE_LENGTH' o wartości 1,000,000, ograniczającą maksymalną długość generowanej sekwencji

    #ORIGINAL
    #length = int(input("Podaj długość sekwencji: "))
    #MODIFIED(Zapewnienie możliwości pracy z błędami powstałymi w wyniku nieprawidłowego wprowadzenia długości sekwencji.)
    # Pobiera długość sekwencji od użytkownika z walidacją
    while True: # Rozpoczyna pętlę nieskończoną, która będzie kontynuowana dopóki użytkownik nie poda poprawnej długości
        try: # Rozpoczyna blok 'try' do obsługi potencjalnych błędów podczas konwersji danych wejściowych
            length = int(input("Podaj długość sekwencji: ")) # Prosi użytkownika o podanie długości sekwencji i próbuje przekonwertować ją na liczbę całkowitą
            if length <= 0: # Sprawdza, czy podana długość jest liczbą dodatnią
                raise ValueError("Długość musi być liczbą dodatnią.") # Jeśli nie, zgłasza błąd typu ValueError z odpowiednim komunikatem
            if length > MAX_SEQUENCE_LENGTH: # Sprawdza, czy podana długość nie przekracza maksymalnej dozwolonej długości
                raise ValueError(f"Długość sekwencji nie może przekroczyć {MAX_SEQUENCE_LENGTH}.") # Jeśli tak, zgłasza błąd typu ValueError
            break # Jeśli długość jest poprawna, przerywa pętlę 'while'
        except ValueError as e: # Przechwytuje błąd typu ValueError (np. gdy użytkownik wpisze tekst zamiast liczby lub liczba jest nieprawidłowa)
            print(f"Błąd: {e}. Spróbuj ponownie.") # Wyświetla komunikat o błędzie i pętla kontynuuje iterację, prosząc ponownie o dane

    # Pobiera dane od użytkownika
    sequence_id = input("Podaj nazwę (ID) sekwencji: ") # Prosi użytkownika o podanie ID sekwencji i zapisuje je w zmiennej 'sequence_id'
    description = input("Podaj opis sekwencji: ") # Prosi użytkownika o podanie opisu sekwencji i zapisuje go w zmiennej 'description'
    name = input("Podaj imię: ") # Prosi użytkownika o podanie imienia, które zostanie wstawione do sekwencji, i zapisuje je w zmiennej 'name'

    # Generuje losową sekwencję DNA
    dna_sequence = generate_dna_sequence(length) # Wywołuje funkcję 'generate_dna_sequence' z podaną przez użytkownika długością, aby stworzyć losową sekwencję DNA

    # Wstawia imię do sekwencji
    sequence_with_name = insert_name_in_sequence(dna_sequence, name) # Wywołuje funkcję 'insert_name_in_sequence', aby wstawić podane imię do oryginalnej sekwencji DNA (choć ta zmienna używana jest tylko do zapisu do FASTA, statystyki są liczone z oryginalnej)

    # Oblicza statystyki (na podstawie oryginalnej sekwencji bez imienia)
    stats = calculate_statistics(dna_sequence) # Wywołuje funkcję 'calculate_statistics' na oryginalnej sekwencji DNA ('dna_sequence') bez wstawionego imienia, aby obliczyć statystyki nukleotydów

    #ORIGINAL
    #for key, value in stats.items():  # Iteruje przez każdą parę klucz-wartość w słowniku 'stats'
    #    if key == 'C/G to A/T ratio':  # Sprawdza, czy klucz to stosunek GC/AT
    #        print(f"%Ratio GC/AT: {value:.2f}")  # Wyświetla stosunek GC/AT na konsoli
    #    else:  # Dla pozostałych statystyk (procentowa zawartość nukleotydów)
    #       print(f"{key}: {value:.2f}%")  # Wyświetla procentową zawartość nukleotydu na konsoli
    #MODIFIED(Dodanie funkcjonalności umożliwiającej zapisywanie wyników obliczeń w pliku w celu ich późniejszego przeglądania)
    # Zapisuje statystyki do pliku oraz wyświetla je
    print("\nStatystyki sekwencji DNA:") # Wyświetla nagłówek dla statystyk
    file_name_statistics = f"{sequence_id}_statistics.txt" # Tworzy nazwę pliku dla statystyk, np. 'ID_sekwencji_statistics.txt'
    try: # Rozpoczyna blok 'try' do obsługi potencjalnych błędów podczas zapisu do pliku
        with open(file_name_statistics, "w") as stats_file: # Otwiera plik do zapisu statystyk
            for key, value in stats.items(): # Iteruje przez każdą parę klucz-wartość w słowniku 'stats'
                if key == 'C/G to A/T ratio': # Sprawdza, czy klucz to stosunek GC/AT
                    stats_file.write(f"Ratio GC/AT: {value:.2f}\n") # Zapisuje stosunek GC/AT do pliku, sformatowany do dwóch miejsc po przecinku
                    print(f"%Ratio GC/AT: {value:.2f}") # Wyświetla stosunek GC/AT na konsoli
                else: # Dla pozostałych statystyk (procentowa zawartość nukleotydów)
                    stats_file.write(f"{key}: {value:.2f}%\n") # Zapisuje procentową zawartość nukleotydu do pliku, sformatowaną do dwóch miejsc po przecinku
                    print(f"{key}: {value:.2f}%") # Wyświetla procentową zawartość nukleotydu na konsoli
        print(f"Statystyki zapisane do pliku {file_name_statistics}.") # Informuje o pomyślnym zapisaniu pliku ze statystykami
    except IOError as e: # Przechwytuje błędy wejścia/wyjścia (np. brak uprawnień do zapisu)
        print(f"Nie udało się zapisać statystyk: {e}") # Wyświetla komunikat o błędzie zapisu statystyk

    # Zapisuje sekwencję do pliku FASTA
    save_to_fasta(sequence_id, description, sequence_with_name) # Wywołuje funkcję 'save_to_fasta', aby zapisać sekwencję DNA (z wstawionym imieniem) do pliku w formacie FASTA

    # Rysuje wykres statystyk
    plot_nucleotide_statistics(stats) # Wywołuje funkcję 'plot_nucleotide_statistics', aby stworzyć i wyświetlić wykres słupkowy statystyk nukleotydów

#ORIGINAL
#brak funkcjonalności
#MODIFIED(Teraz jest możliwość zobaczyć graficzną reprezentację wyników)
def plot_nucleotide_statistics(stats):# Definiuje funkcję programu o nazwie 'plot_nucleotide_statistics' dla rysowania grafiku
    nucleotides = ['A', 'C', 'G', 'T']# Wybieramy tylko klucze odpowiadające nukleotydom (A, C, G, T), ignorujemy 'C/G to A/T ratio'

    values = [stats[n] for n in nucleotides]# Tworzymy listę wartości procentowych dla każdego nukleotydu według kolejności w 'nucleotides'

    plt.figure(figsize=(8, 6))# Tworzymy nową figurę (obszar wykresu) o rozmiarze 8x6 cali
    plt.bar(nucleotides, values, color=['blue', 'green', 'orange', 'red'])# Tworzymy wykres słupkowy (bar chart) z kolorami przypisanymi do nukleotydów
    plt.title('Procentowy udział nukleotydów w sekwencji DNA')# Ustawiamy tytuł wykresu
    plt.xlabel('Nukleotyd')# Etykieta osi X
    plt.ylabel('Zawartość (%)')# Etykieta osi Y
    plt.ylim(0, 100)# Zakres osi Y od 0 do 100 (procenty)
    plt.grid(axis='y', linestyle='--', alpha=0.7)# Dodajemy poziome linie siatki dla lepszej czytelności (linia przerywana, przezroczystość 0.7)

    for i, v in enumerate(values):# Dla każdego słupka dodajemy tekst z wartością procentową nad nim
        plt.text(i, v + 1, f"{v:.2f}%", ha='center', fontweight='bold')# Współrzędne: i — pozycja słupka, v + 1 — trochę nad słupkiem; zaokrąglenie do 2 miejsc

    plt.tight_layout()# Automatyczne dopasowanie elementów wykresu, żeby się nie nakładały
    plt.show()# Wyświetlamy gotowy wykres w oknie graficznym



if __name__ == "__main__": # Standardowa konstrukcja w Pythonie; sprawdza, czy skrypt jest uruchamiany bezpośrednio (a nie importowany jako moduł)
    main()  # Uruchamia główną funkcję programu