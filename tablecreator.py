#!/usr/bin/env python3
"""
Скрипт для создания датасета из данных метилирования GEO
"""

import requests
import re
import time
import random
import pandas as pd
import numpy as np
import argparse
import sys
import os
import xml.etree.ElementTree as ET

class MethylationDatasetCreator:
    def __init__(self, output_dir: str = "methylation_dataset"):
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)
        self.session = requests.Session()
    
    def find_methylation_datasets(self, max_datasets: int = 1000):
        keywords = [
            "DNA methylation",
            "methylation array", 
            "bisulfite sequencing",
            "Illumina Methylation",
            "450k",
            "850k", 
            "EPIC array",
            "CpG methylation",
            "methylome"
        ]

        all_gse_ids = set()

        for keyword in keywords:
            for page in range(1, 21):
                try:
                    params = {
                        'term': f'({keyword}) AND "gse"[Filter]',
                        'retmax': 100,
                        'retstart': (page - 1) * 100,
                        'db': 'gds',
                        'retmode': 'json'
                    }

                    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
                    response = self.session.get(url, params=params, timeout=30)
                    data = response.json()

                    gse_ids = data.get('esearchresult', {}).get('idlist', [])
                    
                    if not gse_ids:
                        break

                    new_ids = [f"GSE{id}" for id in gse_ids]
                    all_gse_ids.update(new_ids)

                    if len(all_gse_ids) >= max_datasets:
                        all_gse_ids = set(list(all_gse_ids)[:max_datasets])
                        break

                    time.sleep(random.uniform(0.3, 1.0))

                except Exception:
                    continue

        return list(all_gse_ids)[:max_datasets]

    def get_dataset_summary(self, gse_id: str):
        try:
            params = {
                'db': 'gds',
                'id': gse_id.replace('GSE', ''),
                'retmode': 'xml'
            }
            
            url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
            response = self.session.get(url, params=params, timeout=30)
            root = ET.fromstring(response.content)
            
            summary_data = {
                'gse_id': gse_id,
                'title': '',
                'summary': '',
                'organism': '',
                'platform': '',
                'samples_count': 0
            }
            
            for item in root.findall('.//Item'):
                if item.get('Name') == 'title':
                    summary_data['title'] = item.text or ''
                elif item.get('Name') == 'summary':
                    summary_data['summary'] = item.text or ''
                elif item.get('Name') == 'taxon':
                    summary_data['organism'] = item.text or ''
                elif item.get('Name') == 'platform':
                    summary_data['platform'] = item.text or ''
                elif item.get('Name') == 'samples':
                    summary_data['samples_count'] = int(item.text) if item.text else 0
            
            return summary_data
            
        except Exception:
            return None

    def extract_sample_metadata(self, gse_id: str):
        try:
            params = {
                'db': 'gds',
                'id': gse_id.replace('GSE', ''),
                'retmode': 'xml'
            }
            
            url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
            response = self.session.get(url, params=params, timeout=30)
            root = ET.fromstring(response.content)
            
            samples_metadata = []
            samples = []
            
            for item in root.findall('.//Item'):
                if item.get('Name') == 'samples':
                    samples = [sample.text for sample in item.findall('Item')]
                    break
            
            for sample_id in samples[:3]:
                try:
                    sample_meta = self._get_sample_details(sample_id)
                    if sample_meta:
                        samples_metadata.append(sample_meta)
                    time.sleep(0.1)
                except Exception:
                    continue
            
            if samples_metadata:
                return pd.DataFrame(samples_metadata)
            else:
                return None
                
        except Exception:
            return None

    def _get_sample_details(self, sample_id: str):
        try:
            params = {
                'db': 'gds',
                'id': sample_id,
                'retmode': 'xml'
            }
            
            url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
            response = self.session.get(url, params=params, timeout=30)
            root = ET.fromstring(response.content)
            
            sample_meta = {
                'sample_id': sample_id,
                'title': '',
                'organism': '',
                'age': np.nan,
                'sex': 'unknown',
                'tissue': 'unknown'
            }
            
            characteristics = []
            for item in root.findall('.//Item'):
                if item.get('Name') == 'title':
                    sample_meta['title'] = item.text or ''
                elif item.get('Name') == 'taxon':
                    sample_meta['organism'] = item.text or ''
                elif item.get('Name') == 'characteristics':
                    characteristics.append(item.text or '')
            
            sample_meta.update(self._parse_characteristics(characteristics))
            return sample_meta
            
        except Exception:
            return None

    def _parse_characteristics(self, characteristics: list):
        demographics = {
            'age': np.nan,
            'sex': 'unknown', 
            'tissue': 'unknown'
        }
        
        for char in characteristics:
            if not char:
                continue
                
            char_lower = char.lower()
            
            age_match = re.search(r'age[:\s]*(\d+\.?\d*)', char_lower)
            if age_match:
                try:
                    demographics['age'] = float(age_match.group(1))
                except (ValueError, TypeError):
                    pass
            
            if 'male' in char_lower or 'gender: m' in char_lower:
                demographics['sex'] = 'male'
            elif 'female' in char_lower or 'gender: f' in char_lower:
                demographics['sex'] = 'female'
            
            tissues = ['blood', 'brain', 'liver', 'skin', 'lung', 'heart']
            for tissue in tissues:
                if tissue in char_lower:
                    demographics['tissue'] = tissue
                    break
        
        return demographics

    def create_dataset(self, max_datasets: int = 1000):
        gse_ids = self.find_methylation_datasets(max_datasets)
        
        if not gse_ids:
            return False
        
        self._save_gse_list(gse_ids)
        
        datasets_info = []
        for gse_id in gse_ids[:100]:
            dataset_summary = self.get_dataset_summary(gse_id)
            if dataset_summary:
                datasets_info.append(dataset_summary)
            time.sleep(0.3)
        
        all_samples_metadata = []
        for gse_id in gse_ids[:20]:
            samples_df = self.extract_sample_metadata(gse_id)
            if samples_df is not None:
                samples_df['gse_id'] = gse_id
                all_samples_metadata.append(samples_df)
            time.sleep(0.3)
        
        self._save_final_datasets(datasets_info, all_samples_metadata)
        return True

    def _save_gse_list(self, gse_ids: list):
        df_ids = pd.DataFrame({'GSE_ID': gse_ids})
        df_ids.to_csv(os.path.join(self.output_dir, 'methylation_gse_ids.csv'), index=False)

    def _save_final_datasets(self, datasets_info: list, all_samples_metadata: list):
        if datasets_info:
            df_datasets = pd.DataFrame(datasets_info)
            df_datasets.to_csv(os.path.join(self.output_dir, 'methylation_datasets_info.csv'), index=False)
        
        if all_samples_metadata:
            df_samples = pd.concat(all_samples_metadata, ignore_index=True)
            df_samples.to_csv(os.path.join(self.output_dir, 'methylation_samples_metadata.csv'), index=False)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--max-datasets', '-m', type=int, default=1000)
    parser.add_argument('--output', '-o', default='methylation_dataset')
    
    args = parser.parse_args()
    
    creator = MethylationDatasetCreator(args.output)
    creator.create_dataset(max_datasets=args.max_datasets)

if __name__ == "__main__":
    main()
