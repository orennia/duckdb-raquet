#!/usr/bin/env python3
"""Convert a binary file to a C array (replacement for xxd -i)."""
import sys
import os

def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <input_file> <output_file>", file=sys.stderr)
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    var_name = os.path.basename(input_file).replace('.', '_')

    with open(input_file, 'rb') as f:
        data = f.read()

    with open(output_file, 'w') as f:
        f.write(f'unsigned char {var_name}[] = {{\n')
        for i in range(0, len(data), 12):
            chunk = data[i:i+12]
            line = ', '.join(f'0x{b:02x}' for b in chunk)
            if i + 12 < len(data):
                line += ','
            f.write(f'  {line}\n')
        f.write('};\n')
        f.write(f'unsigned int {var_name}_len = {len(data)};\n')

    print(f"Generated {output_file} ({len(data)} bytes)")

if __name__ == '__main__':
    main()
