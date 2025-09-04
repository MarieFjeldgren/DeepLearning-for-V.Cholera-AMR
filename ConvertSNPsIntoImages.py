#Imports
import os
import math
import numpy as np
from PIL import Image

# Selecting colors for nucleotides
NUC_COLORS = {
    'A': (0, 128, 0),    # green
    'C': (0, 0, 200),    # blue
    'G': (255, 0, 0),    # red
    'T': (255, 165, 0),  # orange
    'N': (105, 105, 105) # gray
}
PAD_COLOR = (0, 0, 0)  # black padding for blanks


def fasta_to_records(path):
    """Parse FASTA file into a list"""
    records = []
    header, seq = None, []
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    records.append((header, "".join(seq)))
                header = line[1:]
                seq = []
            else:
                seq.append(line)
        if header:
            records.append((header, "".join(seq)))
    return records


def seq_to_image(seq):
    """Convert the SNP nucleotide sequence to a near-square RGB image."""
    seq = seq.upper()
    n = len(seq)
    width = math.ceil(math.sqrt(n))
    height = math.ceil(n / width)
    total = width * height

    img = np.zeros((height, width, 3), dtype=np.uint8)

    for i, base in enumerate(seq):
        r, c = divmod(i, width)
        color = NUC_COLORS.get(base, NUC_COLORS['N'])
        img[r, c] = color

    # Fill padding with black
    for i in range(len(seq), total):
        r, c = divmod(i, width)
        img[r, c] = PAD_COLOR

    return Image.fromarray(img, 'RGB')


# ---- Run section ----
output_folder = "fasta_images"     # <-- all PNGs will be saved here
fasta_file = "snp.aln.fa"  # <-- input FASTA file

# Create folder if it doesn’t exist
os.makedirs(output_folder, exist_ok=True)

records = fasta_to_records(fasta_file)

lengths = [len(seq) for _, seq in records]
if len(set(lengths)) != 1:
    raise ValueError(f"Sequences have different lengths: {set(lengths)}")
else:
    print(f"All sequences have length {lengths[0]}")


for header, seq in records:
    img = seq_to_image(seq)
    safe_name = "".join(c if c.isalnum() else "_" for c in header[:30])  # clean filename
    filename = os.path.join(output_folder, f"{safe_name}.png")
    img.save(filename)






