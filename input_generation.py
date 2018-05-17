import numpy as np
import random
import argparse


def get_argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_file', type=str, default="input.txt")
    parser.add_argument('--num_chunks', type=int, default=100)            
    return parser

def generate_random_data(data_size):
    data_list = [];
    for i in range(data_size):
        alphabet = 'abcdefghijklmnopqrstuvwxyz'
        data_list.append(random.choice(alphabet))
    return data_list


def write_data_to_file(data_list, file_name):
    thefile = open(file_name, 'w');
    for char in data_list:
        thefile.write("%s" %char)

def main():
    parser = get_argument_parser()
    config = parser.parse_args()
    data_list = generate_random_data(32*config.num_chunks)
    write_data_to_file(data_list, config.input_file)

if __name__ == '__main__':
    main()
