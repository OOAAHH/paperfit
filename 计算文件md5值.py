import sys
import hashlib

def calculate_hash(file_path):
    # Create a new MD5 hash object
    hash_obj = hashlib.md5()
    with open(file_path, 'rb') as file:
        # Read the file in 4096 byte chunks
        for chunk in iter(lambda: file.read(4096), b''):
            hash_obj.update(chunk)
    # Return the hexadecimal representation of the hash
    return hash_obj.hexdigest()

def write_hash_to_file(input_file_path, hash_value):
    names = input_file_path
    names = names.split('/')[-1]
    outfiles = f'{names}_hashes_8T.txt'
    with open(outfiles, 'a') as file:
        file.write(f"{input_file_path}:{hash_value}\n")

'''if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py /path/file")
        sys.exit(1)
    file_path = sys.argv[1]
    hash_value = calculate_hash(file_path)
    write_hash_to_file("hashes.txt", hash_value)
    print(f"Hash value for {file_path} written to hashes.txt")'''

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py /path/file")
        sys.exit(1)
    input_file_path = sys.argv[1]
    hash_value = calculate_hash(input_file_path)
    write_hash_to_file(input_file_path, hash_value)
    print(f"Hash value for {input_file_path} written to hashes.txt")
