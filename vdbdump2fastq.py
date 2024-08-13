import re
import logging
import gzip

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def parse_linkage_group(linkage_group):
    """Extract CB and UB from linkage group"""
    try:
        cb_match = re.search(r'CB:(\w+)', linkage_group)
        ub_match = re.search(r'UB:(\w+)', linkage_group)
        cb = cb_match.group(1) if cb_match else None
        ub = ub_match.group(1) if ub_match else None
        return cb, ub
    except Exception as e:
        logging.error(f"Error parsing linkage group: {e}")
        return None, None

def quality_to_phred33(quality_str):
    """Convert quality numbers to Phred+33 ASCII characters"""
    try:
        qualities = map(int, quality_str.split(", "))
        return ''.join(chr(q + 33) for q in qualities)
    except Exception as e:
        logging.error(f"Error converting quality scores: {e}")
        return None

def process_lines(lines):
    """Process lines to extract data into a group"""
    current_group = {}
    fastq_entries = []
    total_groups = 0
    processed_count = 0

    for line in lines:
        line = line.strip()

        # Skip empty lines
        if not line:
            continue

        if line.startswith("ALIGNMENT_COUNT:"):
            if current_group:
                total_groups += 1
                fastq_entry = generate_fastq_entry(current_group)
                if fastq_entry:
                    fastq_entries.append(fastq_entry)
                    processed_count += 1
            current_group = {}

        key_value = line.split(":", 1)
        if len(key_value) == 2:
            key = key_value[0].strip()
            value = key_value[1].strip()
            current_group[key] = value

    # Process the last group if it exists
    if current_group:
        total_groups += 1
        fastq_entry = generate_fastq_entry(current_group)
        if fastq_entry:
            fastq_entries.append(fastq_entry)
            processed_count += 1

    logging.info(f"Total groups found: {total_groups}")
    logging.info(f"Total groups processed: {processed_count}")
    logging.info(f"Total valid reads generated: {len(fastq_entries)}")

    return fastq_entries

def generate_fastq_entry(group):
    """Generate a FASTQ entry from a data group"""
    try:
        name = group.get('NAME')
        read = group.get('READ')
        quality = group.get('QUALITY')
        linkage_group = group.get('LINKAGE_GROUP')
        spot_group = group.get('SPOT_GROUP')

        if name and read and quality and linkage_group and spot_group:
            quality = quality_to_phred33(quality)
            cb, ub = parse_linkage_group(linkage_group)
            if cb and ub:
                rg = spot_group
                fastq_name = f"@{name}:{rg}_{cb}_{ub}"
                return f"{fastq_name}\n{read}\n+\n{quality}\n"
        return None
    except Exception as e:
        logging.error(f"Error generating FASTQ entry: {e}")
        return None

def write_fastq_file(fastq_entries, output_file):
    """Write the fastq entries to a file"""
    try:
        with open(output_file, 'w') as f:
            f.writelines(fastq_entries)
        logging.info(f"FASTQ file '{output_file}' written successfully.")
    except Exception as e:
        logging.error(f"Error writing to file: {e}")

def read_input_file(input_file):
    """Read the input file, supporting both plain text and gzip files"""
    try:
        if input_file.endswith('.gz'):
            with gzip.open(input_file, 'rt') as f:
                lines = f.readlines()
        else:
            with open(input_file, 'r') as f:
                lines = f.readlines()
        return lines
    except Exception as e:
        logging.error(f"Error reading input file: {e}")
        return []

if __name__ == "__main__":
    import sys
    input_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else "output.fastq"

    # Read the input file, handling .gz files if necessary
    lines = read_input_file(input_file)

    # Process lines to extract fastq entries
    fastq_entries = process_lines(lines)

    # Write to fastq file
    write_fastq_file(fastq_entries, output_file)
