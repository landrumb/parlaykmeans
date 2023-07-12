import csv
import struct

def bin_to_csv(bin_file, csv_file):
    with open(bin_file, 'rb') as binfile, open(csv_file, 'w', newline='') as csvfile:
        header_size = struct.calcsize('i')
        header = struct.unpack('i', binfile.read(header_size))[0]
        data = binfile.read()

        csv_data = []
        pos = 0
        while pos < len(data):
            row_size = struct.unpack_from('i', data, pos)[0]
            pos += struct.calcsize('i')
            row_data = []
            for _ in range(row_size):
                item_size = struct.unpack_from('i', data, pos)[0]
                pos += struct.calcsize('i')
                item = struct.unpack_from(f'{item_size}s', data, pos)[0].decode()
                pos += item_size
                row_data.append(item)
            csv_data.append(row_data)

        writer = csv.writer(csvfile)
        writer.writerows(csv_data)

bin_file_path = input("Enter the BIN file path: ")
csv_file_path = input("Enter the CSV file path: ")

bin_to_csv(bin_file_path, csv_file_path)
