#!/usr/bin/env python3

from ancseq.params import Params


pm = Params('ancseq')
args = pm.set_options()

import os
import sys
import gzip
import shutil
import subprocess as sbp
from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import SimpleFastaParser
from ancseq.utils import time_stamp, clean_cmd, call_log


class ancseq(object):

    def __init__(self, args):
        self.args = args

    def built_tree(self):
        print(time_stamp(),
              'IQ-TREE로 계통수 구축 중...',
              flush=True)
        out_dir_00 = os.path.join(self.args.out, '00_tree')
        seq = os.path.join(out_dir_00, os.path.basename(self.args.seq))
        os.mkdir(out_dir_00)
        shutil.copyfile(self.args.seq, seq)
        cmd = f'iqtree -s {seq} \
                       -st {self.args.mode} \
                       -T {self.args.threads} \
                       -m {self.args.model}'
        if self.args.fast:
            cmd += f' --fast \
                      --alrt {self.args.bootstrap}'
        else:
            cmd += f' -B {self.args.bootstrap}'
        if self.args.outgroup != None:
            cmd += f' -o {self.args.outgroup}'
        cmd += f' 1> {out_dir_00}/00_iqtree.out \
                  2> {out_dir_00}/00_iqtree.err'
        cmd = clean_cmd(cmd)
        try:
            sbp.run(cmd,
                    stdout=sbp.DEVNULL,
                    stderr=sbp.DEVNULL,
                    shell=True,
                    check=True)
        except sbp.CalledProcessError:
            call_log(f'{out_dir_00}/00_iqtree.err', cmd)
            sys.exit(1)
        print(time_stamp(),
              'IQ-TREE 작업이 성공적으로 완료되었습니다.',
              flush=True)
        
    def check_best_model(self):
        iqtree_log = os.path.join(self.args.out, '00_tree', f'{os.path.basename(self.args.seq)}.log')
        with open(iqtree_log) as log:
            for line in log:
                if line.startswith('Best-fit model:'): # 원본 로그 파일 형식 유지
                    self.args.model = line.split()[2]
                    print(time_stamp(),
                        f'IQ-TREE에서 최적 모델은 {self.args.model} 입니다.',
                        flush=True)
    
    def reconstruct_ancestral_state(self):
        print(time_stamp(),
              '조상 상태 재구축 중...',
              flush=True)
        out_dir_10 = os.path.join(self.args.out, '10_asr')
        seq = os.path.join(out_dir_10, os.path.basename(self.args.seq))
        tree = os.path.join(self.args.out, '00_tree', f'{os.path.basename(self.args.seq)}.treefile')
        os.mkdir(out_dir_10)
        shutil.copyfile(self.args.seq, seq)
        cmd = f'iqtree -asr \
                       -s {seq} \
                       -te {tree} \
                       -st {self.args.mode} \
                       -T {self.args.threads} \
                       -m {self.args.model} \
                       -keep_empty_seq'
        if self.args.outgroup != None:
            cmd += f' -o {self.args.outgroup}'
        cmd += f' 1> {out_dir_10}/10_iqtree.out \
                  2> {out_dir_10}/10_iqtree.err'
        cmd = clean_cmd(cmd)
        try:
            sbp.run(cmd,
                    stdout=sbp.DEVNULL,
                    stderr=sbp.DEVNULL,
                    shell=True,
                    check=True)
        except sbp.CalledProcessError:
            call_log(f'{out_dir_10}/10_iqtree.err', cmd)
            sys.exit(1)        
        print(time_stamp(),
              '조상 서열 재구축이 성공적으로 완료되었습니다.',
              flush=True)
    
    def alignment_to_binary(self):
        self.binary_file_name = os.path.join(self.args.out, '20_indels', f'{os.path.basename(self.args.seq)}.binary')
        dna_missing_chars = {'N', 'X', '?', '.', 'O', '~', '!'}  # 집합으로 변경하여 검색 성능 향상
        aa_missing_chars = {'X', '?', '.', '~', '!'} # 집합으로 변경하여 검색 성능 향상

        with open(self.args.seq) as fasta_file, open(self.binary_file_name, 'w') as binary_file:
            for header, seq in SimpleFastaParser(fasta_file):
                bin_seq = []  # 문자열 이어붙이기 대신 리스트 사용 후 join
                if self.args.mode == 'CODON':
                    for i in range(0, len(seq), 3):
                        codon = seq[i:i+3]
                        if codon == '---':
                            bin_seq.append('1')
                        elif any(c in dna_missing_chars for c in codon): # any와 생성기 표현식 사용
                            bin_seq.append('-')
                        else:
                            bin_seq.append('0')
                elif self.args.mode == 'DNA':
                    for char_seq in seq:
                        if char_seq == '-':
                            bin_seq.append('1')
                        elif char_seq in dna_missing_chars:
                            bin_seq.append('-')
                        else:
                            bin_seq.append('0')
                elif self.args.mode == 'AA':
                    for char_seq in seq:
                        if char_seq == '-':
                            bin_seq.append('1')
                        elif char_seq in aa_missing_chars:
                            bin_seq.append('-')
                        else:
                            bin_seq.append('0')
                else:
                    print(time_stamp(),
                          f'잘못된 모드 {self.args.mode}가 주어졌습니다.',
                          file=sys.stderr,
                          flush=True)
                    sys.exit(1)
                binary_file.write(f'>{header}\n')
                binary_file.write(''.join(bin_seq) + '\n')
        
    def reconstruct_indels(self):
        print(time_stamp(),
              'INDEL 재구축 중...',
              flush=True)
        out_dir_20 = os.path.join(self.args.out, '20_indels')
        tree = os.path.join(self.args.out, '00_tree', f'{os.path.basename(self.args.seq)}.treefile')
        os.mkdir(out_dir_20)
        self.alignment_to_binary()
        cmd = f'iqtree -asr \
                       -s {self.binary_file_name} \
                       -te {tree} \
                       -st BIN \
                       -T {self.args.threads} \
                       -blfix \
                       -keep_empty_seq \
                       -m JC2'
        if self.args.outgroup != None:
            cmd += f' -o {self.args.outgroup}'
        cmd += f' 1> {out_dir_20}/20_iqtree.out \
                  2> {out_dir_20}/20_iqtree.err'
        cmd = clean_cmd(cmd)
        try:
            sbp.run(cmd,
                    stdout=sbp.DEVNULL,
                    stderr=sbp.DEVNULL,
                    shell=True,
                    check=True)
        except sbp.CalledProcessError:
            call_log(f'{out_dir_20}/20_iqtree.err', cmd)
            sys.exit(1)        
        print(time_stamp(),
              'INDEL 재구축이 성공적으로 완료되었습니다.',
              flush=True)

    def merge_results(self):
        out_dir_30 = os.path.join(self.args.out, '30_result')
        os.makedirs(out_dir_30, exist_ok=True)

        path_probs_result_10 = os.path.join(self.args.out, '10_asr', f'{os.path.basename(self.args.seq)}.state')
        path_probs_result_20 = os.path.join(self.args.out, '20_indels', f'{os.path.basename(self.args.seq)}.binary.state')
        path_merged = f'{out_dir_30}/ancestral_state_result.tsv'
        path_fasta_with_gap = f'{out_dir_30}/ancestral_state_result_with_gap.fasta'
        path_fasta = f'{out_dir_30}/ancestral_state_result.fasta'

        gap_char = '---' if self.args.mode == 'CODON' else '-'
        previous_node = None
        node_seq_list = []

        try:
            with open(path_probs_result_10) as f_probs10, \
                 open(path_probs_result_20) as f_probs20, \
                 open(path_merged, 'w') as f_merged, \
                 open(path_fasta_with_gap, 'w') as f_fasta_gap, \
                 open(path_fasta, 'w') as f_fasta:

                for line_10, line_20 in zip(f_probs10, f_probs20):
                    if line_10.startswith('#'):
                        continue

                    cols_10 = line_10.rstrip('\n').split('\t')
                    if line_10.startswith('Node\tSite'):
                        f_merged.write('Node\tSite\tState\t{}\t{}\n'.format('\t'.join([p.replace('p_', '') for p in cols_10[3:]]), gap_char))
                        previous_node = None
                        node_seq_list = []
                        continue

                    cols_20 = line_20.rstrip('\n').split('\t')
                    node_10, site_10, state = cols_10[0], cols_10[1], cols_10[2]
                    node_20, site_20 = cols_20[0], cols_20[1]

                    if (node_10, site_10) != (node_20, site_20):
                        print(time_stamp(), f'{path_probs_result_10} 와 {path_probs_result_20}의 노드/사이트 ({node_10}, {site_10}) vs ({node_20}, {site_20})가 일치하지 않습니다.', file=sys.stderr, flush=True)
                        sys.exit(1)

                    if previous_node != node_10 and previous_node is not None:
                        f_fasta_gap.write(f'>{previous_node}\n{"".join(node_seq_list)}\n')
                        f_fasta.write(f'>{previous_node}\n{"".join(node_seq_list).replace("-", "")}\n')
                        node_seq_list = []

                    probs_str = [format(round(float(p), 4), '.5f') for p in cols_10[3:]]
                    p_1 = float(cols_20[4])

                    current_state_char = gap_char if p_1 > self.args.min_gap_prob else state
                    f_merged.write('{}\t{}\t{}\t{}\t{:.4f}\n'.format(node_10, site_10, current_state_char, '\t'.join(probs_str), p_1))
                    node_seq_list.append(current_state_char)
                    previous_node = node_10

                if previous_node is not None and node_seq_list:
                    f_fasta_gap.write(f'>{previous_node}\n{"".join(node_seq_list)}\n')
                    f_fasta.write(f'>{previous_node}\n{"".join(node_seq_list).replace("-", "")}\n')

        except FileNotFoundError as e:
            print(time_stamp(), f"오류: 파일 없음 - {e.filename}", file=sys.stderr, flush=True)
            sys.exit(1)
        except Exception as e:
            print(time_stamp(), f"결과 병합 중 오류 발생: {e}", file=sys.stderr, flush=True)
            sys.exit(1)


    def calculate_codon_prob(self):
        print(time_stamp(), '코돈 확률 계산 중...', flush=True)
        out_dir_30 = os.path.join(self.args.out, '30_result')
        merged_file_path = f'{out_dir_30}/ancestral_state_result.tsv'
        codon_prob_file_path = f'{out_dir_30}/ancestral_state_result.codon_prob.tsv'

        from itertools import product # itertools import 추가
        nucleotides = ['A', 'C', 'G', 'T']
        all_codons = [''.join(c) for c in product(nucleotides, repeat=3)] # product 사용

        try:
            with open(merged_file_path) as merged_file, open(codon_prob_file_path, 'w') as codon_prob_file:
                header = merged_file.readline().rstrip('\n').split('\t')
                # Correctly extract nucleotide probability headers (e.g., A, C, G, T)
                # This assumes the order in merged_file is always A, C, G, T, followed by gap prob
                nuc_prob_headers = header[3:7]


                codon_prob_file.write(f"Node\tSite\tState\t{'\t'.join(all_codons)}\t---\n")

                nuc_dicts = [{}, {}, {}] # Store nucleotide probabilities for a codon
                gap_probs = [0.0, 0.0, 0.0] # Store gap probabilities for a codon

                for line in merged_file:
                    line = line.rstrip('\n')
                    cols = line.split('\t')

                    node, site_str, state_from_merged = cols[0], cols[1], cols[2]
                    site = int(site_str)
                    pos_in_codon = (site -1) % 3 # 0, 1, or 2

                    # Store probabilities for the current nucleotide position
                    # Ensure columns 3-7 contain A, C, G, T probabilities and column 7 is gap
                    current_nuc_probs = {nuc: float(p) for nuc, p in zip(nuc_prob_headers, cols[3:7])}
                    nuc_dicts[pos_in_codon] = current_nuc_probs
                    gap_probs[pos_in_codon] = float(cols[7])

                    if pos_in_codon == 2: # End of a codon
                        codon_probs_calculated = {
                            codon: nuc_dicts[0].get(codon[0],0) * nuc_dicts[1].get(codon[1],0) * nuc_dicts[2].get(codon[2],0)
                            for codon in all_codons
                        }
                        # Ensure all_codons are present, default to 0 if a nucleotide was missing (e.g. from a gap in previous position)
                        for c in all_codons:
                            if c not in codon_probs_calculated:
                                codon_probs_calculated[c] = 0.0

                        codon_probs_calculated['---'] = sum(gap_probs) / 3.0

                        prob_values_str_list = []
                        for c in all_codons: # Use defined order
                            prob_values_str_list.append(format(round(codon_probs_calculated.get(c, 0.0), 5), '.5f')) # Use .get for safety
                        prob_values_str_list.append(format(round(codon_probs_calculated.get('---',0.0), 5), '.5f'))


                        sorted_codon_probs = sorted(
                            [(codon, prob) for codon, prob in codon_probs_calculated.items() if codon != '---'],
                            key=lambda x: x[1],
                            reverse=True
                        )

                        display_state = '---' if codon_probs_calculated['---'] >= self.args.min_gap_prob else (sorted_codon_probs[0][0] if sorted_codon_probs else 'NNN')

                        codon_prob_file.write(f"{node}\t{int(site/3)+1}\t{display_state}\t{'\t'.join(prob_values_str_list)}\n")
                        nuc_dicts = [{}, {}, {}]
                        gap_probs = [0.0, 0.0, 0.0]

        except FileNotFoundError as e:
            print(time_stamp(), f"오류: 파일 없음 - {e.filename}", file=sys.stderr, flush=True)
            sys.exit(1)
        except Exception as e:
            print(time_stamp(), f"코돈 확률 계산 중 오류 발생: {e} (Line: {sys.exc_info()[-1].tb_lineno})", file=sys.stderr, flush=True)
            sys.exit(1)

        print(time_stamp(), '코돈 확률 계산이 완료되었습니다.', flush=True)

    def _get_top_states(self, probs_dict, max_report, min_report_prob, is_codon_mode, min_gap_threshold):
        """Helper function to get top states (codons or AAs) based on probabilities."""
        gap_char = '---' if is_codon_mode else '-'

        if probs_dict.get(gap_char, 0.0) >= min_gap_threshold:
            top_states = [gap_char] + (max_report - 1) * ['']
            top_state_probs_str = [format(round(probs_dict.get(gap_char, 0.0), 5), '.5f')] + \
                                  (max_report - 1) * ['']
            if is_codon_mode:
                top_aas = ['-'] + (max_report - 1) * ['']
                top_aa_probs_str = [format(round(probs_dict.get(gap_char, 0.0), 5), '.5f')] + \
                                   (max_report - 1) * ['']
        else:
            # Remove gap for sorting if it exists and is below threshold
            filtered_probs = {k: v for k, v in probs_dict.items() if k != gap_char}

            sorted_states = sorted(filtered_probs.items(), key=lambda x: x[1], reverse=True)

            top_states = [k if v >= min_report_prob else '' for k, v in sorted_states][:max_report]
            top_states.extend([''] * (max_report - len(top_states))) # Ensure list has max_report elements

            top_state_probs_str = [format(round(v, 5), '.5f') if v >= min_report_prob else '' for _, v in sorted_states][:max_report]
            top_state_probs_str.extend([''] * (max_report - len(top_state_probs_str)))

            if is_codon_mode:
                aa_probs = {}
                for codon, prob in filtered_probs.items(): # Use filtered_probs
                    try:
                        # Ensure Seq object is created correctly for translation
                        aa = str(Seq(codon).translate(to_stop=True)) # Use to_stop=True to handle stop codons as '*'
                        if aa == '*': aa = '-' # Or however you want to represent stop codons as AA
                    except Exception:
                        aa = 'X' # Represent unknown translation as 'X'
                    aa_probs[aa] = aa_probs.get(aa, 0.0) + prob

                sorted_aas = sorted(aa_probs.items(), key=lambda x: x[1], reverse=True)
                top_aas = [k if v >= min_report_prob else '' for k, v in sorted_aas][:max_report]
                top_aas.extend([''] * (max_report - len(top_aas)))

                top_aa_probs_str = [format(round(v, 5), '.5f') if v >= min_report_prob else '' for _, v in sorted_aas][:max_report]
                top_aa_probs_str.extend([''] * (max_report - len(top_aa_probs_str)))

        return (top_states, top_state_probs_str, top_aas, top_aa_probs_str) if is_codon_mode else (top_states, top_state_probs_str)


    def sort_ambiguous_states(self):
        """Sorts ambiguous codons or amino acids based on mode and writes to a TSV file."""
        is_codon_mode = self.args.mode in ['DNA', 'CODON']
        mode_specific_message = '코돈' if is_codon_mode else '아미노산'
        print(time_stamp(), f"모호한 {mode_specific_message} 정렬 중...", flush=True)

        out_dir_30 = os.path.join(self.args.out, '30_result')
        
        if self.args.mode == 'DNA':
            input_file_path = f"{out_dir_30}/ancestral_state_result.codon_prob.tsv"
        elif self.args.mode == 'CODON' or self.args.mode == 'AA':
             input_file_path = f"{out_dir_30}/ancestral_state_result.tsv"
        else:
            print(time_stamp(), f"sort_ambiguous_states: 지원되지 않는 모드 {self.args.mode}", file=sys.stderr, flush=True)
            return

        output_file_path = f"{out_dir_30}/ancestral_state_result.sort.tsv"

        try:
            with open(input_file_path) as infile, open(output_file_path, 'w') as outfile:
                header_line = infile.readline().rstrip('\n')
                header_cols = header_line.split('\t')

                # Determine state keys from header (e.g. AAA, AAC ... or A, C, D ...)
                # For codon_prob.tsv, states are from col 3 up to (but not including) the last '---' col
                # For ancestral_state_result.tsv (CODON/AA), states are from col 3 up to the last col (gap prob)
                if self.args.mode == 'DNA': # Reading from ancestral_state_result.codon_prob.tsv
                    # Header: Node Site State <codon1> <codon2> ... <codon64> ---
                    state_keys = header_cols[3:-1] # All codons
                else: # Reading from ancestral_state_result.tsv (for CODON or AA)
                    # Header: Node Site State <state1> <state2> ... <stateN> <gap_char>
                    state_keys = header_cols[3:-1] # All AAs or Codons

                gap_key_in_header = header_cols[-1] # This should be '---' or '-'

                for line in infile:
                    cols = line.rstrip('\n').split('\t')
                    node, site_val_str = cols[0], cols[1]

                    # Create dictionary of state -> probability
                    current_probs = {key: float(val) for key, val in zip(state_keys, cols[3:3+len(state_keys)])}
                    # Add gap probability
                    current_probs[gap_key_in_header] = float(cols[-1])


                    if is_codon_mode: # CODON or DNA mode
                        top_codons, top_codon_probs_str, top_aas, top_aa_probs_str = self._get_top_states(
                            current_probs, self.args.max_report, self.args.min_prob, True, self.args.min_gap_prob
                        )
                        outfile.write(f"{node}\t{site_val_str}\t{'\t'.join(top_codons)}\t{'\t'.join(top_codon_probs_str)}\t{'\t'.join(top_aas)}\t{'\t'.join(top_aa_probs_str)}\n")
                    else: # AA mode
                        top_aas_list, top_aa_probs_str_list = self._get_top_states(
                            current_probs, self.args.max_report, self.args.min_prob, False, self.args.min_gap_prob
                        )
                        outfile.write(f"{node}\t{site_val_str}\t{'\t'.join(top_aas_list)}\t{'\t'.join(top_aa_probs_str_list)}\n")

        except FileNotFoundError:
            print(time_stamp(), f"오류: 입력 파일 없음 - {input_file_path}", file=sys.stderr, flush=True)
            sys.exit(1)
        except Exception as e:
            print(time_stamp(), f"모호한 상태 정렬 중 오류 발생: {e} (Line: {sys.exc_info()[-1].tb_lineno})", file=sys.stderr, flush=True)

            sys.exit(1)

        print(time_stamp(), f"모호한 {mode_specific_message}이(가) 정렬되었습니다.", flush=True)

    def copy_treefile(self):
        out_dir_10 = os.path.join(self.args.out, '10_asr')
        out_dir_30 = os.path.join(self.args.out, '30_result')
        treefile = os.path.join(out_dir_10, f'{os.path.basename(self.args.seq)}.treefile')
        shutil.copyfile(treefile, f'{out_dir_30}/ancestral_state_result.treefile')
        
    def gzip_table(self):
        out_dir_10 = os.path.join(self.args.out, '10_asr')
        out_dir_20 = os.path.join(self.args.out, '20_indels')
        out_dir_30 = os.path.join(self.args.out, '30_result')
        statefile_10 = os.path.join(out_dir_10, f'{os.path.basename(self.args.seq)}.state')
        statefile_20 = os.path.join(out_dir_20, f'{os.path.basename(self.args.seq)}.binary.state')
        marged_prefix = os.path.join(out_dir_30, 'ancestral_state_result')
        if os.path.isfile(statefile_10):
            with open(statefile_10, 'rb') as f_in:
                with gzip.open(f'{statefile_10}.gz', 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            os.remove(statefile_10)
        if os.path.isfile(statefile_20):
            with open(statefile_20, 'rb') as f_in:
                with gzip.open(f'{statefile_20}.gz', 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            os.remove(statefile_20)
        if os.path.isfile(f'{marged_prefix}.tsv'):
            with open(f'{marged_prefix}.tsv', 'rb') as f_in:
                with gzip.open(f'{marged_prefix}.tsv.gz', 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            os.remove(f'{marged_prefix}.tsv')
        if os.path.isfile(f'{marged_prefix}.codon_prob.tsv'):
            with open(f'{marged_prefix}.codon_prob.tsv', 'rb') as f_in:
                with gzip.open(f'{marged_prefix}.codon_prob.tsv.gz', 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            os.remove(f'{marged_prefix}.codon_prob.tsv')

    def run(self):
        print(time_stamp(), 'ancseq 실행 시작.', flush=True)
        os.mkdir(self.args.out)
        if self.args.asr_only:
            print(time_stamp(),
                  '계통수 구축을 건너뛰었습니다.',
                  flush=True)
            out_dir_00 = os.path.join(self.args.out, '00_tree')
            # Ensure 00_tree directory exists
            os.makedirs(out_dir_00, exist_ok=True)
            seq_target_path = os.path.join(out_dir_00, os.path.basename(self.args.seq))
            shutil.copyfile(self.args.seq, seq_target_path)

            # Check if source treefile exists before copying
            source_treefile = f'{self.args.seq}.treefile'
            if os.path.isfile(source_treefile):
                shutil.copyfile(source_treefile, f'{seq_target_path}.treefile')
            else:
                print(time_stamp(), f"경고: --asr-only가 지정되었지만, 원본 트리 파일({source_treefile})을 찾을 수 없습니다.", file=sys.stderr, flush=True)

            source_logfile = f'{self.args.seq}.log'
            if os.path.isfile(source_logfile): # Check if source logfile exists
                shutil.copyfile(source_logfile, f'{seq_target_path}.log')
        else:
            self.built_tree()

        # Check best model only if tree was built or log file was copied
        potential_log_file = os.path.join(self.args.out, '00_tree', f'{os.path.basename(self.args.seq)}.log')
        if self.args.model == 'MFP' and os.path.isfile(potential_log_file):
            self.check_best_model()

        self.reconstruct_ancestral_state()
        self.reconstruct_indels()
        self.merge_results()

        if self.args.mode == 'DNA' and not self.args.stop_codon_prob:
            self.calculate_codon_prob()
            self.sort_ambiguous_states() # 통합된 함수 호출
        elif self.args.mode == 'CODON':
            self.sort_ambiguous_states() # 통합된 함수 호출
        elif self.args.mode == 'AA':
            self.sort_ambiguous_states() # 통합된 함수 호출

        self.copy_treefile()
        self.gzip_table()
        print(time_stamp(), 'ancseq가 성공적으로 완료되었습니다!', flush=True)


def main():
    # pm and args are already defined at the top level
    # No, they should be inside main or passed if ancseq is a reusable class outside this script context.
    # For this script, current setup is fine as it's the main entry point.
    ancseq_runner = ancseq(args)
    ancseq_runner.run()

if __name__ == '__main__':
    main()



