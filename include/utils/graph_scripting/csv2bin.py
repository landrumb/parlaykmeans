import csv
import struct

def csv_to_bin(csv_file, bin_file):
    with open(csv_file, 'r') as csvfile, open(bin_file, 'wb') as binfile:
        reader = csv.reader(csvfile)
        rows = list(reader)
        binfile.write(struct.pack('i', len(rows)))
        for row in rows:
            binfile.write(struct.pack('i', len(row)))
            for item in row:
                item_bytes = item.encode()
                binfile.write(struct.pack('i', len(item_bytes)))
                binfile.write(item_bytes)

csv_file_path = input("Enter the CSV file path: ")
bin_file_path = input("Enter the BIN file path: ")

csv_to_bin(csv_file_path, bin_file_path)
