#!/usr/bin/env python3

import argparse
import csv
import json
import os


def parse_kma_results(kma_results_path):
    kma_results = []

    int_fields = [
        'score',
        'expected',
        'template_length',
    ]

    float_fields = [
        'template_identity',
        'template_coverage',
        'query_identity',
        'query_coverage',
        'depth',
        'q_value',
        'p_value',
    ]
    
    with open(kma_results_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            for field in int_fields:
                row[field] = int(row[field])
            for field in float_fields:
                row[field] = float(row[field])
            kma_results.append(row)
            
    return kma_results


def find_seq_in_db(db_path, kma_template_name):
    matches = []
    with open(db_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                defline = line.strip()[1:]
                defline_without_whitespace = defline.replace(' ', '')
                if kma_template_name == defline_without_whitespace:
                    matches.append(defline)

    return matches
                
            


def main(args):
    kma_results = parse_kma_results(args.kma_results)
    kma_results_sorted_by_q_value_descending = list(sorted(kma_results, key=lambda x: x['q_value'], reverse=True))
    if len(kma_results_sorted_by_q_value_descending) > 0:
        best_match_kma_template = kma_results_sorted_by_q_value_descending[0]['template']
        matches = find_seq_in_db(args.db, best_match_kma_template)
    
        for match in matches:
            print(match)
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('kma_results')
    parser.add_argument('--db')
    args = parser.parse_args()
    main(args)
