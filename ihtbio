import time
from Bio import Entrez
from Bio import SeqIO

# --- НАСТРОЙКА ENTREZ ---
Entrez.email = 'ваша_почта@gmail.com' # Обязательно укажите свою почту
time.sleep(1 / 10) # Защита от блокировок (ограничение частоты запросов)

print("--- Задание 1: Загрузка генома SARS-CoV-2 ---")
# Используем корректный ID NC_045512.2 (актуальная версия в NCBI)
handle = Entrez.efetch(db='nucleotide', id='NC_045512.2', rettype='gb', retmode='text')
record = SeqIO.read(handle, 'genbank')
handle.close()

seq = str(record.seq).upper()
print(f"Геном успешно загружен. Длина последовательности: {len(seq)} нуклеотидов.\n")


print("--- Задание 2: Поиск потенциальных ORF областей ---")
# Пороги длины в аминокислотах (зададим, например, от 100 до 1000 аа)
min_length = 100 
max_length = 1000
STOP_CODONS = {'TAA', 'TAG', 'TGA'}

orfs = []

for frame in range(3):
    i = frame
    while i < len(seq) - 2:
        if seq[i:i+3] == 'ATG':
            start = i
            j = i + 3
            while j < len(seq) - 2:
                codon = seq[j:j+3]
                if codon in STOP_CODONS:
                    break
                j += 3

            if j < len(seq) - 2 and seq[j:j+3] in STOP_CODONS:
                length_aa = (j - start) // 3
                if max_length >= length_aa >= min_length:
                    orfs.append({
                        "frame": frame + 1,
                        "start": start + 1,
                        "end": j + 3,
                        "length_aa": length_aa
                    })
            i += 3
        else:
            i += 3

print(f"Найдено потенциальных ORF (длиной {min_length}-{max_length} aa): {len(orfs)}")
for orf in orfs[:5]: # Выведем первые 5 для примера
    print(orf)
print("...\n")


print("--- Задание 3: Вывод аннотированных CDS из GenBank ---")
for feature in record.features:
    if feature.type == 'CDS':
        name = feature.qualifiers.get('gene', ['?'])[0]
        product = feature.qualifiers.get('product', ['?'])[0]
        start = int(feature.location.start)
        end = int(feature.location.end)
        cds_seq = record.seq[start:end]
        print(f"Ген {name} ({product}): координаты {start}–{end}, длина {len(cds_seq)} нт")
