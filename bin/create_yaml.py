import csv
import sys
import yaml

def csv_to_yaml(csv_path, output_path, output_dir):
    data = {}
    with open(csv_path, 'r') as csv_file:
        reader = csv.reader(csv_file)
        next(reader)  # Skip header row
        for row in reader:
            sample_name = row[0]
            condition = row[1]
            replicate = row[2]
            data_path = sample_name
            
            if condition not in data:
                data[condition] = {}
            
            data[condition][replicate] = data_path
    
    output = {'data': data, 'out': output_dir}
    
    with open(output_path, 'w') as yaml_file:
        yaml.dump(output, yaml_file, default_flow_style=False)

# Example usage
csv_file_path = sys.argv[1]
yaml_output_path = 'xpore_config.yaml'
output_directory = '.'
csv_to_yaml(csv_file_path, yaml_output_path, output_directory)