# how to test for a different protein:
# - download xsl file from CRISPOR (https://crispor.gi.ucsc.edu/) for choosen cds sequence and convert it to csv
# - add CRISPOR.csv to bio_test_files
# - add csv result file downloaded from this website to bio_test_files
# use def test_ABE_P04439(self): and change the folowing lines:
# - df_crispor = pd.read_csv('designer/tests/bio_test_files/ YOUR CRISPOR FILE.csv')
#   df_base_editing = pd.read_csv('designer/tests/bio_test_files/ YOUR BASE EDITOR FILE.csv')
# - change Exon boundaries, cds sequence and ensembl transcript ID for your protein
# - test with killed_by_only_ABE/ killed_by_only_CBE if testing only ABE/ CBE else use none
# - test with killed_by_exon_filter
# - test with killed_by_no_AS_mutation_ABE/ killed_by_no_AS_mutation_CBE if testing only ABE/CBE else both

from django.test import TestCase
import unittest
from unittest.mock import patch, MagicMock
from django.test import TestCase
import pandas as pd
import math
import re
from Bio.Seq import Seq


class BioTest(TestCase):

    @staticmethod
    def print_data(crispor_seqs, base_editor_seqs):
        nur_in_crispor = crispor_seqs - base_editor_seqs

        print(f"CRISPOR hat insgesamt: {len(crispor_seqs)} Guides")
        print(f"Unser Tool hat insgesamt: {len(base_editor_seqs)} Guides")
        print(f"\nGuides die CRISPOR hat, wir NICHT: {len(nur_in_crispor)}")
        print(f"Guides die wir haben, CRISPOR NICHT: {len(base_editor_seqs - crispor_seqs)}")
        print(f"Guides die beide haben: {len(crispor_seqs & base_editor_seqs)}")

    @staticmethod
    def parse_crispor_position(guide_id):
        """Extrahiert Position und Strand aus CRISPOR guideId z.B. '624rev' oder '305forw'"""
        match = re.match(r'(\d+)(rev|forw)', str(guide_id))
        if match:
            return int(match.group(1)), match.group(2)
        return None, None

    @staticmethod
    def find_pos_in_base_editor(cds, pos, sgRNAseq, strand):

        pam = sgRNAseq[-3:]

        if strand == 'rev':
            pam = str(Seq(pam).complement())
            pam = ''.join(reversed(pam))

        pam_start = 0
        pam_end = 0


        for offset in [-1, 0, 1]:
            start = pos + offset
            if cds[start:start+3] == pam:
                pam_start = start
                pam_end = start+3
                break
        else:
            raise ValueError(f"Found no PAM at {pos} /n {cds[pos-1:pos+2]}, {cds[pos:pos+3]}, {cds[pos+1:pos+4]} for {pam}")

        editing_window = sgRNAseq[3:8] # 4 - 8
        positions = []

        # 0 <- PAM
        if(strand == 'forw'):
            editing_window_start = pam_end - 23 + 3
            editing_window_end = pam_end - 23 + 8

            if(cds[editing_window_start: editing_window_end] != editing_window):
                raise ValueError(f"Found wrong Editing Window: {cds[editing_window_start: editing_window_end]}, {editing_window} for {pos}")

            for i in range(editing_window_start+1, editing_window_end+3):
                if((i-3) % 3 == 0):
                    positions.append(int((i)/3))

        else: # strand == 'rev'

            editing_window = str(Seq(editing_window).complement())
            editing_window = ''.join(reversed(editing_window))

            editing_window_start = pam_end + 23 - 8 - 3
            editing_window_end = pam_end + 23 - 3 - 3

            if(cds[editing_window_start: editing_window_end] != editing_window):
                raise ValueError(f"Found wrong Editing Window: {cds[editing_window_start: editing_window_end]}({editing_window_start},{editing_window_end}){pam_end}, {editing_window} at {pos}")

            for i in range(editing_window_start+1, editing_window_end+3):
                if((i-3) % 3 == 0):
                    positions.append(int(i/3))

        return (positions,(editing_window_start,editing_window_end))

    @staticmethod
    def killed_by_only_ABE(row):

        seq = row['targetSeq'].upper()
        window =  seq[3:8]

        if 'A' in window:
            return False
        else:
            return True

    @staticmethod
    def killed_by_only_CBE(row):

        seq = row['targetSeq'].upper()
        window =  seq[3:8]

        if 'C' in window:
            return False
        else:
            return True

    @staticmethod
    def killed_by_exon_filter(row, exon_boundaries, cds):
        near_boundary = False

        pam = row['targetSeq'][-3:]
        pos = row['pos']

        if row['strand'] == 'rev':
            pam = str(Seq(pam).complement())
            pam = ''.join(reversed(pam))

        pam_start = 0
        pam_end = 0


        for offset in [-1, 0, 1]:
            start = pos + offset
            if cds[start:start+3] == pam:
                pam_start = start
                pam_end = start+3
                break
        else:
            raise ValueError(f"Found no PAM at {pos} /n {cds[pos-1:pos+2]}, {cds[pos:pos+3]}, {cds[pos+1:pos+4]} for {pam}")

        for boundary in exon_boundaries:
            if row['strand'] == 'rev':
                if boundary - 23 <= pam_end <= boundary:
                    near_boundary = True

            elif row['strand'] == 'forw':
                if boundary <= pam_end <= boundary + 23:
                    near_boundary = True

        return near_boundary

    @staticmethod
    def killed_by_no_AS_mutation_ABE(row, position_tuple, cds):

        positions = position_tuple[0]
        window = position_tuple[1]

        strand = row['strand']

        for i in positions:
            codon_start = (i*3 - 3)
            codon = cds[codon_start:codon_start+3]

            # safety check
            if len(codon) != 3:
                continue

            aa = str(Seq(codon).translate())

            mutated_codon = codon
            new_aa = ''

            if(codon_start == window[0] - 2):
                if(strand == 'forw' and codon[2] == 'A'):
                    mutated_codon = codon[:-1] + 'G'
                elif(strand == 'rev' and codon[2] == 'T'):
                    mutated_codon = codon[:-1] + 'C'

            elif(codon_start == window[0] - 1):
                s_list = list(mutated_codon)

                for i in range(max(0, len(s_list) - 2), len(s_list)):
                    if s_list[i] == 'A' and strand == 'forw':
                        s_list[i] = 'G'
                    if strand == 'rev' and s_list[i] == 'T':
                        s_list[i] = 'C'

                mutated_codon = ''.join(s_list)

            elif(codon_start == window[1] - 2):
                s_list = list(mutated_codon)

                for i in range(max(0, len(s_list) - 2), len(s_list)):
                    if(s_list[i] == 'A' and strand == 'forw'):
                        s_list[i] = 'G'
                    if(strand == 'rev' and s_list[i] == 'T'):
                        s_list[i] = 'C'

                mutated_codon = ''.join(s_list)

            elif(codon_start == window[1] - 1):
                if(strand == 'forw' and codon[0] == 'A'):
                    mutated_codon = 'G' + codon[1:]
                if(strand == 'rev' and codon[0] == 'T'):
                    mutated_codon = 'C' + codon[1:]

            else: # whole codon in editing window
                if(strand == 'forw'):

                    mutated_codon = codon.replace('A', 'G')

                if(strand == 'rev'):
                    mutated_codon = codon.replace('T', 'C')

            new_aa = str(Seq(mutated_codon).translate())

            stop = str(Seq('TAA').translate())

            if new_aa != aa and new_aa != stop:
                return False

        # if we never found a non-silent mutation
        return True  # killed

    @staticmethod
    def killed_by_no_AS_mutation_CBE(row, position_tuple, cds):

        positions = position_tuple[0]
        window = position_tuple[1]

        strand = row['strand']

        for i in positions:
            codon_start = (i*3 - 3)
            codon = cds[codon_start:codon_start+3]

            # safety check
            if len(codon) != 3:
                continue

            aa = str(Seq(codon).translate())

            mutated_codon = codon
            new_aa = ''

            if(codon_start == window[0] - 2):
                if(strand == 'forw' and codon[2] == 'C'):
                    mutated_codon = codon[:-1] + 'T'
                elif(strand == 'rev' and codon[2] == 'G'):
                    mutated_codon = codon[:-1] + 'A'

            elif(codon_start == window[0] - 1):
                s_list = list(mutated_codon)

                for i in range(max(0, len(s_list) - 2), len(s_list)):
                    if s_list[i] == 'C' and strand == 'forw':
                        s_list[i] = 'T'
                    if strand == 'rev' and s_list[i] == 'G':
                        s_list[i] = 'A'

                mutated_codon = ''.join(s_list)

            elif(codon_start == window[1] - 2):
                s_list = list(mutated_codon)

                for i in range(max(0, len(s_list) - 2), len(s_list)):
                    if(s_list[i] == 'C' and strand == 'forw'):
                        s_list[i] = 'T'
                    if(strand == 'rev' and s_list[i] == 'G'):
                        s_list[i] = 'A'

                mutated_codon = ''.join(s_list)

            elif(codon_start == window[1] - 1):
                if(strand == 'forw' and codon[0] == 'C'):
                    mutated_codon = 'T' + codon[1:]
                if(strand == 'rev' and codon[0] == 'G'):
                    mutated_codon = 'A' + codon[1:]

            else: # whole codon in editing window
                if(strand == 'forw'):
                    mutated_codon = codon.replace('C', 'T')

                if(strand == 'rev'):
                    mutated_codon = codon.replace('G', 'A')

            new_aa = str(Seq(mutated_codon).translate())

        # if we never found a non-silent mutation
        return True  # killed

    @staticmethod
    def killed_by_no_AS_mutation_CBE(row, position_tuple, cds):
        if(row['pos'] == 215):
            raise ValueError (f"{position_tuple}")
        positions = position_tuple[0]
        window = position_tuple[1]

        strand = row['strand']

        for i in positions:
            codon_start = (i*3 - 3)
            codon = cds[codon_start:codon_start+3]
            
            # safety check
            if len(codon) != 3:
                continue

            aa = str(Seq(codon).translate())

            mutated_codon = codon
            new_aa = ""

            if(codon_start == window[0] - 2):
                if(strand == 'forw' and codon[2] == 'C'):
                    mutated_codon = codon[:-1] + "T"
                elif(strand == 'rev' and codon[2] == 'G'):
                    mutated_codon = codon[:-1] + "A"

            elif(codon_start == window[0] - 1):
                s_list = list(mutated_codon)

                for i in range(max(0, len(s_list) - 2), len(s_list)):
                    if s_list[i] == 'C' and strand == 'forw':
                        s_list[i] = 'T'
                    if strand == 'rev' and s_list[i] == 'G':
                        s_list[i] = 'A'

                mutated_codon = ''.join(s_list)

            elif(codon_start == window[1] - 2):
                s_list = list(mutated_codon)
               
                for i in range(max(0, len(s_list) - 2), len(s_list)):
                    if(s_list[i] == 'C' and strand == 'forw'):
                        s_list[i] = 'T'
                    if(strand == 'rev' and s_list[i] == 'G'):
                        s_list[i] = 'A'
                
                mutated_codon = ''.join(s_list)
                
            elif(codon_start == window[1] - 1):
                if(strand == 'forw' and codon[0] == 'C'):
                    mutated_codon = "T" + codon[1:]
                if(strand == 'rev' and codon[0] == 'G'):
                    mutated_codon = "A" + codon[1:]
                
            else: # codon ist komplett im editing window
                if(strand == 'forw'):
                    mutated_codon = codon.replace("C", "T")
                    
                if(strand == 'rev'):
                    mutated_codon = codon.replace("G", "A")
            
            new_aa = str(Seq(mutated_codon).translate())  

            stop = str(Seq('TAA').translate())

            if new_aa != aa and new_aa != stop:
                return False

        # if we never found a non-silent mutation
        return True  # killed

# ===========================
# Protein tests
# ===========================

    def test_ABE_P04439(self):
        # test details: P04439, ABE, editing window 4-8, show doubles
        from designer.pipeline import get_exon_boundaries
        from designer.services import run_analysis

        df_crispor = pd.read_csv('designer/tests/bio_test_files/crispor_results_P04439.csv')
        df_base_editing = pd.read_csv('designer/tests/bio_test_files/our_base_results_P04439_ABE.csv')

        ensembl_transcript_ID = 'ENST00000376809'
        cds = 'ATGGCCGTCATGGCGCCCCGAACCCTCCTCCTGCTACTCTCGGGGGCCCTGGCCCTGACCCAGACCTGGGCGGGCTCCCACTCCATGAGGTATTTCTTCACATCCGTGTCCCGGCCCGGCCGCGGGGAGCCCCGCTTCATCGCCGTGGGCTACGTGGACGACACGCAGTTCGTGCGGTTCGACAGCGACGCCGCGAGCCAGAGGATGGAGCCGCGGGCGCCGTGGATAGAGCAGGAGGGGCCGGAGTATTGGGACCAGGAGACACGGAATGTGAAGGCCCAGTCACAGACTGACCGAGTGGACCTGGGGACCCTGCGCGGCTACTACAACCAGAGCGAGGCCGGTTCTCACACCATCCAGATAATGTATGGCTGCGACGTGGGGTCGGACGGGCGCTTCCTCCGCGGGTACCGGCAGGACGCCTACGACGGCAAGGATTACATCGCCCTGAACGAGGACCTGCGCTCTTGGACCGCGGCGGACATGGCGGCTCAGATCACCAAGCGCAAGTGGGAGGCGGCCCATGAGGCGGAGCAGTTGAGAGCCTACCTGGATGGCACGTGCGTGGAGTGGCTCCGCAGATACCTGGAGAACGGGAAGGAGACGCTGCAGCGCACGGACCCCCCCAAGACACATATGACCCACCACCCCATCTCTGACCATGAGGCCACCCTGAGGTGCTGGGCCCTGGGCTTCTACCCTGCGGAGATCACACTGACCTGGCAGCGGGATGGGGAGGACCAGACCCAGGACACGGAGCTCGTGGAGACCAGGCCTGCAGGGGATGGAACCTTCCAGAAGTGGGCGGCTGTGGTGGTGCCTTCTGGAGAGGAGCAGAGATACACCTGCCATGTGCAGCATGAGGGTCTGCCCAAGCCCCTCACCCTGAGATGGGAGCTGTCTTCCCAGCCCACCATCCCCATCGTGGGCATCATTGCTGGCCTGGTTCTCCTTGGAGCTGTGATCACTGGAGCTGTGGTCGCTGCCGTGATGTGGAGGAGGAAGAGCTCAGATAGAAAAGGAGGGAGTTACACTCAGGCTGCAAGCAGTGACAGTGCCCAGGGCTCTGATGTGTCCCTCACAGCTTGTAAAGTGTGA'
        exon_boundaries = [73, 343, 619, 895, 1012, 1045, 1093]


        crispor_seqs = set(df_crispor['targetSeq'].str.upper().str.strip())
        base_editor_seqs = set(df_base_editing['Guide RNA'].str.upper().str.strip())

        self.print_data(crispor_seqs, base_editor_seqs)

        only_in_crispor =  crispor_seqs - base_editor_seqs
        only_in_base_editor = base_editor_seqs - crispor_seqs

        if len(only_in_base_editor) > 0:
            self.fail('Base editor produced sgRNAs not found in CRISPOR')

        df_not_in_base_editor = df_crispor[df_crispor['targetSeq'].str.upper().str.strip().isin(only_in_crispor)].copy()

        df_not_in_base_editor['pos'] = df_not_in_base_editor['#guideId'].apply(lambda x: self.parse_crispor_position(x)[0])
        df_not_in_base_editor['strand'] = df_not_in_base_editor['#guideId'].apply(lambda x: self.parse_crispor_position(x)[1])

        unexpected = []

        for _, row in df_not_in_base_editor.iterrows():
            pos_in_base_editor = self.find_pos_in_base_editor(cds, row['pos'], row['targetSeq'], row['strand'])

            if self.killed_by_only_ABE(row):
                continue
            elif self.killed_by_exon_filter(row, exon_boundaries, cds):
                continue
            elif self.killed_by_no_AS_mutation_ABE(row, pos_in_base_editor, cds):
                continue
            else:
                unexpected.append(row)

        self.assertEqual(len(unexpected), 0, f"Unexplained guides: {len(unexpected)}{unexpected}")

    def test_CBE_P04439(self):
        # test details: P04439, CBE, editing window 4-8, show doubles
        from designer.pipeline import get_exon_boundaries
        from designer.services import run_analysis

        df_crispor = pd.read_csv('designer/tests/bio_test_files/crispor_results_P04439.csv')
        df_base_editing = pd.read_csv('designer/tests/bio_test_files/our_base_results_P04439_CBE.csv')


        ensembl_transcript_ID = 'ENST00000376809'
        cds = 'ATGGCCGTCATGGCGCCCCGAACCCTCCTCCTGCTACTCTCGGGGGCCCTGGCCCTGACCCAGACCTGGGCGGGCTCCCACTCCATGAGGTATTTCTTCACATCCGTGTCCCGGCCCGGCCGCGGGGAGCCCCGCTTCATCGCCGTGGGCTACGTGGACGACACGCAGTTCGTGCGGTTCGACAGCGACGCCGCGAGCCAGAGGATGGAGCCGCGGGCGCCGTGGATAGAGCAGGAGGGGCCGGAGTATTGGGACCAGGAGACACGGAATGTGAAGGCCCAGTCACAGACTGACCGAGTGGACCTGGGGACCCTGCGCGGCTACTACAACCAGAGCGAGGCCGGTTCTCACACCATCCAGATAATGTATGGCTGCGACGTGGGGTCGGACGGGCGCTTCCTCCGCGGGTACCGGCAGGACGCCTACGACGGCAAGGATTACATCGCCCTGAACGAGGACCTGCGCTCTTGGACCGCGGCGGACATGGCGGCTCAGATCACCAAGCGCAAGTGGGAGGCGGCCCATGAGGCGGAGCAGTTGAGAGCCTACCTGGATGGCACGTGCGTGGAGTGGCTCCGCAGATACCTGGAGAACGGGAAGGAGACGCTGCAGCGCACGGACCCCCCCAAGACACATATGACCCACCACCCCATCTCTGACCATGAGGCCACCCTGAGGTGCTGGGCCCTGGGCTTCTACCCTGCGGAGATCACACTGACCTGGCAGCGGGATGGGGAGGACCAGACCCAGGACACGGAGCTCGTGGAGACCAGGCCTGCAGGGGATGGAACCTTCCAGAAGTGGGCGGCTGTGGTGGTGCCTTCTGGAGAGGAGCAGAGATACACCTGCCATGTGCAGCATGAGGGTCTGCCCAAGCCCCTCACCCTGAGATGGGAGCTGTCTTCCCAGCCCACCATCCCCATCGTGGGCATCATTGCTGGCCTGGTTCTCCTTGGAGCTGTGATCACTGGAGCTGTGGTCGCTGCCGTGATGTGGAGGAGGAAGAGCTCAGATAGAAAAGGAGGGAGTTACACTCAGGCTGCAAGCAGTGACAGTGCCCAGGGCTCTGATGTGTCCCTCACAGCTTGTAAAGTGTGA'
        exon_boundaries = [73, 343, 619, 895, 1012, 1045, 1093]

        crispor_seqs = set(df_crispor['targetSeq'].str.upper().str.strip())
        base_editor_seqs = set(df_base_editing['Guide RNA'].str.upper().str.strip())

        self.print_data(crispor_seqs, base_editor_seqs)

        only_in_crispor =  crispor_seqs - base_editor_seqs
        only_in_base_editor = base_editor_seqs - crispor_seqs

        if len(only_in_base_editor) > 0:
            self.fail('Base editor produced sgRNAs not found in CRISPOR')

        df_not_in_base_editor = df_crispor[df_crispor['targetSeq'].str.upper().str.strip().isin(only_in_crispor)].copy()

        df_not_in_base_editor['pos'] = df_not_in_base_editor['#guideId'].apply(lambda x: self.parse_crispor_position(x)[0])
        df_not_in_base_editor['strand'] = df_not_in_base_editor['#guideId'].apply(lambda x: self.parse_crispor_position(x)[1])

        unexpected = []

        for _, row in df_not_in_base_editor.iterrows():
            pos_in_base_editor = self.find_pos_in_base_editor(cds, row['pos'], row['targetSeq'], row['strand'])

            if self.killed_by_only_CBE(row):
                continue
            elif self.killed_by_exon_filter(row, exon_boundaries, cds):
                continue
            elif self.killed_by_no_AS_mutation_CBE(row, pos_in_base_editor, cds):
                continue
            else:
                unexpected.append(row)

        self.assertEqual(len(unexpected), 0, f"Unexplained guides: {len(unexpected)}{unexpected}")

    def test_ABE_and_CBE_P04439(self):
        # test details: P04439, ABE&CBE, editing window 4-8, show doubles
        from designer.pipeline import get_exon_boundaries
        from designer.services import run_analysis

        df_crispor = pd.read_csv('designer/tests/bio_test_files/crispor_results_P04439.csv')
        df_base_editing = pd.read_csv('designer/tests/bio_test_files/our_base_results_P04439_ABE_and_CBE.csv')


        ensembl_transcript_ID = 'ENST00000376809'
        cds = 'ATGGCCGTCATGGCGCCCCGAACCCTCCTCCTGCTACTCTCGGGGGCCCTGGCCCTGACCCAGACCTGGGCGGGCTCCCACTCCATGAGGTATTTCTTCACATCCGTGTCCCGGCCCGGCCGCGGGGAGCCCCGCTTCATCGCCGTGGGCTACGTGGACGACACGCAGTTCGTGCGGTTCGACAGCGACGCCGCGAGCCAGAGGATGGAGCCGCGGGCGCCGTGGATAGAGCAGGAGGGGCCGGAGTATTGGGACCAGGAGACACGGAATGTGAAGGCCCAGTCACAGACTGACCGAGTGGACCTGGGGACCCTGCGCGGCTACTACAACCAGAGCGAGGCCGGTTCTCACACCATCCAGATAATGTATGGCTGCGACGTGGGGTCGGACGGGCGCTTCCTCCGCGGGTACCGGCAGGACGCCTACGACGGCAAGGATTACATCGCCCTGAACGAGGACCTGCGCTCTTGGACCGCGGCGGACATGGCGGCTCAGATCACCAAGCGCAAGTGGGAGGCGGCCCATGAGGCGGAGCAGTTGAGAGCCTACCTGGATGGCACGTGCGTGGAGTGGCTCCGCAGATACCTGGAGAACGGGAAGGAGACGCTGCAGCGCACGGACCCCCCCAAGACACATATGACCCACCACCCCATCTCTGACCATGAGGCCACCCTGAGGTGCTGGGCCCTGGGCTTCTACCCTGCGGAGATCACACTGACCTGGCAGCGGGATGGGGAGGACCAGACCCAGGACACGGAGCTCGTGGAGACCAGGCCTGCAGGGGATGGAACCTTCCAGAAGTGGGCGGCTGTGGTGGTGCCTTCTGGAGAGGAGCAGAGATACACCTGCCATGTGCAGCATGAGGGTCTGCCCAAGCCCCTCACCCTGAGATGGGAGCTGTCTTCCCAGCCCACCATCCCCATCGTGGGCATCATTGCTGGCCTGGTTCTCCTTGGAGCTGTGATCACTGGAGCTGTGGTCGCTGCCGTGATGTGGAGGAGGAAGAGCTCAGATAGAAAAGGAGGGAGTTACACTCAGGCTGCAAGCAGTGACAGTGCCCAGGGCTCTGATGTGTCCCTCACAGCTTGTAAAGTGTGA'
        exon_boundaries = [73, 343, 619, 895, 1012, 1045, 1093]


        crispor_seqs = set(df_crispor['targetSeq'].str.upper().str.strip())
        base_editor_seqs = set(df_base_editing['Guide RNA'].str.upper().str.strip())

        self.print_data(crispor_seqs, base_editor_seqs)

        only_in_crispor =  crispor_seqs - base_editor_seqs
        only_in_base_editor = base_editor_seqs - crispor_seqs

        if len(only_in_base_editor) > 0:
            self.fail('Base editor produced sgRNAs not found in CRISPOR')

        df_not_in_base_editor = df_crispor[df_crispor['targetSeq'].str.upper().str.strip().isin(only_in_crispor)].copy()

        df_not_in_base_editor['pos'] = df_not_in_base_editor['#guideId'].apply(lambda x: self.parse_crispor_position(x)[0])
        df_not_in_base_editor['strand'] = df_not_in_base_editor['#guideId'].apply(lambda x: self.parse_crispor_position(x)[1])

        unexpected = []

        for _, row in df_not_in_base_editor.iterrows():
            pos_in_base_editor = self.find_pos_in_base_editor(cds, row['pos'], row['targetSeq'], row['strand'])

            if self.killed_by_exon_filter(row, exon_boundaries, cds):
                continue
            elif self.killed_by_no_AS_mutation_ABE(row, pos_in_base_editor, cds):
                continue
            elif self.killed_by_no_AS_mutation_CBE(row, pos_in_base_editor, cds):
                continue
            else:
                unexpected.append(row)

        self.assertEqual(len(unexpected), 0, f"Unexplained guides: {len(unexpected)}{unexpected}")

    def test_ABE_P31785(self):
        # test details: P31785, ABE, editing window 4-8, show doubles
        from designer.pipeline import get_exon_boundaries
        from designer.services import run_analysis

        df_crispor = pd.read_csv('designer/tests/bio_test_files/crispor_results_P31785.csv')
        df_base_editing = pd.read_csv('designer/tests/bio_test_files/our_base_results_P31785_ABE.csv')

        ensembl_transcript_ID = 'ENST00000374202'
        cds = 'ATGTTGAAGCCATCATTACCATTCACATCCCTCTTATTCCTGCAGCTGCCCCTGCTGGGAGTGGGGCTGAACACGACAATTCTGACGCCCAATGGGAATGAAGACACCACAGCTGATTTCTTCCTGACCACTATGCCCACTGACTCCCTCAGTGTTTCCACTCTGCCCCTCCCAGAGGTTCAGTGTTTTGTGTTCAATGTCGAGTACATGAATTGCACTTGGAACAGCAGCTCTGAGCCCCAGCCTACCAACCTCACTCTGCATTATTGGTACAAGAACTCGGATAATGATAAAGTCCAGAAGTGCAGCCACTATCTATTCTCTGAAGAAATCACTTCTGGCTGTCAGTTGCAAAAAAAGGAGATCCACCTCTACCAAACATTTGTTGTTCAGCTCCAGGACCCACGGGAACCCAGGAGACAGGCCACACAGATGCTAAAACTGCAGAATCTGGTGATCCCCTGGGCTCCAGAGAACCTAACACTTCACAAACTGAGTGAATCCCAGCTAGAACTGAACTGGAACAACAGATTCTTGAACCACTGTTTGGAGCACTTGGTGCAGTACCGGACTGACTGGGACCACAGCTGGACTGAACAATCAGTGGATTATAGACATAAGTTCTCCTTGCCTAGTGTGGATGGGCAGAAACGCTACACGTTTCGTGTTCGGAGCCGCTTTAACCCACTCTGTGGAAGTGCTCAGCATTGGAGTGAATGGAGCCACCCAATCCACTGGGGGAGCAATACTTCAAAAGAGAATCCTTTCCTGTTTGCATTGGAAGCCGTGGTTATCTCTGTTGGCTCCATGGGATTGATTATCAGCCTTCTCTGTGTGTATTTCTGGCTGGAACGGACGATGCCCCGAATTCCCACCCTGAAGAACCTAGAGGATCTTGTTACTGAATACCACGGGAACTTTTCGGCCTGGAGTGGTGTGTCTAAGGGACTGGCTGAGAGTCTGCAGCCAGACTACAGTGAACGACTCTGCCTCGTCAGTGAGATTCCCCCAAAAGGAGGGGCCCTTGGGGAGGGGCCTGGGGCCTCCCCATGCAACCAGCATAGCCCCTACTGGGCCCCCCCATGTTACACCCTAAAGCCTGAAACCTGA'
        exon_boundaries = [115, 269, 454, 594, 757, 854, 924]

        crispor_seqs = set(df_crispor['targetSeq'].str.upper().str.strip())
        base_editor_seqs = set(df_base_editing['Guide RNA'].str.upper().str.strip())

        self.print_data(crispor_seqs, base_editor_seqs)

        only_in_crispor =  crispor_seqs - base_editor_seqs
        only_in_base_editor = base_editor_seqs - crispor_seqs

        if len(only_in_base_editor) > 0:
            self.fail('Base editor produced sgRNAs not found in CRISPOR')

        df_not_in_base_editor = df_crispor[df_crispor['targetSeq'].str.upper().str.strip().isin(only_in_crispor)].copy()

        df_not_in_base_editor['pos'] = df_not_in_base_editor['#guideId'].apply(lambda x: self.parse_crispor_position(x)[0])
        df_not_in_base_editor['strand'] = df_not_in_base_editor['#guideId'].apply(lambda x: self.parse_crispor_position(x)[1])

        unexpected = []

        for _, row in df_not_in_base_editor.iterrows():
            pos_in_base_editor = self.find_pos_in_base_editor(cds, row['pos'], row['targetSeq'], row['strand'])

            if self.killed_by_only_ABE(row):
                continue
            elif self.killed_by_exon_filter(row, exon_boundaries, cds):
                continue
            elif self.killed_by_no_AS_mutation_ABE(row, pos_in_base_editor, cds):
                continue
            else:
                unexpected.append(row)

        self.assertEqual(len(unexpected), 0, f"Unexplained guides: {len(unexpected)}{unexpected}")

    def test_CBE_P31785(self):
        # test details: P31785, CBE, editing window 4-8, show doubles
        from designer.pipeline import get_exon_boundaries
        from designer.services import run_analysis

        df_crispor = pd.read_csv('designer/tests/bio_test_files/crispor_results_P31785.csv')
        df_base_editing = pd.read_csv('designer/tests/bio_test_files/our_base_results_P31785_CBE.csv')

        ensembl_transcript_ID = 'ENST00000374202'
        cds = 'ATGTTGAAGCCATCATTACCATTCACATCCCTCTTATTCCTGCAGCTGCCCCTGCTGGGAGTGGGGCTGAACACGACAATTCTGACGCCCAATGGGAATGAAGACACCACAGCTGATTTCTTCCTGACCACTATGCCCACTGACTCCCTCAGTGTTTCCACTCTGCCCCTCCCAGAGGTTCAGTGTTTTGTGTTCAATGTCGAGTACATGAATTGCACTTGGAACAGCAGCTCTGAGCCCCAGCCTACCAACCTCACTCTGCATTATTGGTACAAGAACTCGGATAATGATAAAGTCCAGAAGTGCAGCCACTATCTATTCTCTGAAGAAATCACTTCTGGCTGTCAGTTGCAAAAAAAGGAGATCCACCTCTACCAAACATTTGTTGTTCAGCTCCAGGACCCACGGGAACCCAGGAGACAGGCCACACAGATGCTAAAACTGCAGAATCTGGTGATCCCCTGGGCTCCAGAGAACCTAACACTTCACAAACTGAGTGAATCCCAGCTAGAACTGAACTGGAACAACAGATTCTTGAACCACTGTTTGGAGCACTTGGTGCAGTACCGGACTGACTGGGACCACAGCTGGACTGAACAATCAGTGGATTATAGACATAAGTTCTCCTTGCCTAGTGTGGATGGGCAGAAACGCTACACGTTTCGTGTTCGGAGCCGCTTTAACCCACTCTGTGGAAGTGCTCAGCATTGGAGTGAATGGAGCCACCCAATCCACTGGGGGAGCAATACTTCAAAAGAGAATCCTTTCCTGTTTGCATTGGAAGCCGTGGTTATCTCTGTTGGCTCCATGGGATTGATTATCAGCCTTCTCTGTGTGTATTTCTGGCTGGAACGGACGATGCCCCGAATTCCCACCCTGAAGAACCTAGAGGATCTTGTTACTGAATACCACGGGAACTTTTCGGCCTGGAGTGGTGTGTCTAAGGGACTGGCTGAGAGTCTGCAGCCAGACTACAGTGAACGACTCTGCCTCGTCAGTGAGATTCCCCCAAAAGGAGGGGCCCTTGGGGAGGGGCCTGGGGCCTCCCCATGCAACCAGCATAGCCCCTACTGGGCCCCCCCATGTTACACCCTAAAGCCTGAAACCTGA'
        exon_boundaries = [115, 269, 454, 594, 757, 854, 924]

        crispor_seqs = set(df_crispor['targetSeq'].str.upper().str.strip())
        base_editor_seqs = set(df_base_editing['Guide RNA'].str.upper().str.strip())

        self.print_data(crispor_seqs, base_editor_seqs)

        only_in_crispor =  crispor_seqs - base_editor_seqs
        only_in_base_editor = base_editor_seqs - crispor_seqs

        if len(only_in_base_editor) > 0:
            self.fail('Base editor produced sgRNAs not found in CRISPOR')

        df_not_in_base_editor = df_crispor[df_crispor['targetSeq'].str.upper().str.strip().isin(only_in_crispor)].copy()

        df_not_in_base_editor['pos'] = df_not_in_base_editor['#guideId'].apply(lambda x: self.parse_crispor_position(x)[0])
        df_not_in_base_editor['strand'] = df_not_in_base_editor['#guideId'].apply(lambda x: self.parse_crispor_position(x)[1])

        unexpected = []

        for _, row in df_not_in_base_editor.iterrows():
            pos_in_base_editor = self.find_pos_in_base_editor(cds, row['pos'], row['targetSeq'], row['strand'])

            if self.killed_by_only_CBE(row):
                continue
            elif self.killed_by_exon_filter(row, exon_boundaries, cds):
                continue
            elif self.killed_by_no_AS_mutation_CBE(row, pos_in_base_editor, cds):
                continue
            else:
                unexpected.append(row)

        self.assertEqual(len(unexpected), 0, f"Unexplained guides: {len(unexpected)}{unexpected}")

    def test_ABE_and_CBE_P31785(self):
        # test details: P31785, ABE&CBE, editing window 4-8, show doubles
        from designer.pipeline import get_exon_boundaries
        from designer.services import run_analysis

        df_crispor = pd.read_csv('designer/tests/bio_test_files/crispor_results_P31785.csv')
        df_base_editing = pd.read_csv('designer/tests/bio_test_files/our_base_results_P31785_ABE_and_CBE.csv')

        ensembl_transcript_ID = 'ENST00000374202'
        cds = 'ATGTTGAAGCCATCATTACCATTCACATCCCTCTTATTCCTGCAGCTGCCCCTGCTGGGAGTGGGGCTGAACACGACAATTCTGACGCCCAATGGGAATGAAGACACCACAGCTGATTTCTTCCTGACCACTATGCCCACTGACTCCCTCAGTGTTTCCACTCTGCCCCTCCCAGAGGTTCAGTGTTTTGTGTTCAATGTCGAGTACATGAATTGCACTTGGAACAGCAGCTCTGAGCCCCAGCCTACCAACCTCACTCTGCATTATTGGTACAAGAACTCGGATAATGATAAAGTCCAGAAGTGCAGCCACTATCTATTCTCTGAAGAAATCACTTCTGGCTGTCAGTTGCAAAAAAAGGAGATCCACCTCTACCAAACATTTGTTGTTCAGCTCCAGGACCCACGGGAACCCAGGAGACAGGCCACACAGATGCTAAAACTGCAGAATCTGGTGATCCCCTGGGCTCCAGAGAACCTAACACTTCACAAACTGAGTGAATCCCAGCTAGAACTGAACTGGAACAACAGATTCTTGAACCACTGTTTGGAGCACTTGGTGCAGTACCGGACTGACTGGGACCACAGCTGGACTGAACAATCAGTGGATTATAGACATAAGTTCTCCTTGCCTAGTGTGGATGGGCAGAAACGCTACACGTTTCGTGTTCGGAGCCGCTTTAACCCACTCTGTGGAAGTGCTCAGCATTGGAGTGAATGGAGCCACCCAATCCACTGGGGGAGCAATACTTCAAAAGAGAATCCTTTCCTGTTTGCATTGGAAGCCGTGGTTATCTCTGTTGGCTCCATGGGATTGATTATCAGCCTTCTCTGTGTGTATTTCTGGCTGGAACGGACGATGCCCCGAATTCCCACCCTGAAGAACCTAGAGGATCTTGTTACTGAATACCACGGGAACTTTTCGGCCTGGAGTGGTGTGTCTAAGGGACTGGCTGAGAGTCTGCAGCCAGACTACAGTGAACGACTCTGCCTCGTCAGTGAGATTCCCCCAAAAGGAGGGGCCCTTGGGGAGGGGCCTGGGGCCTCCCCATGCAACCAGCATAGCCCCTACTGGGCCCCCCCATGTTACACCCTAAAGCCTGAAACCTGA'
        exon_boundaries = [115, 269, 454, 594, 757, 854, 924]

        crispor_seqs = set(df_crispor['targetSeq'].str.upper().str.strip())
        base_editor_seqs = set(df_base_editing['Guide RNA'].str.upper().str.strip())

        self.print_data(crispor_seqs, base_editor_seqs)

        only_in_crispor =  crispor_seqs - base_editor_seqs
        only_in_base_editor = base_editor_seqs - crispor_seqs

        if len(only_in_base_editor) > 0:
            self.fail('Base editor produced sgRNAs not found in CRISPOR')

        df_not_in_base_editor = df_crispor[df_crispor['targetSeq'].str.upper().str.strip().isin(only_in_crispor)].copy()

        df_not_in_base_editor['pos'] = df_not_in_base_editor['#guideId'].apply(lambda x: self.parse_crispor_position(x)[0])
        df_not_in_base_editor['strand'] = df_not_in_base_editor['#guideId'].apply(lambda x: self.parse_crispor_position(x)[1])

        unexpected = []

        for _, row in df_not_in_base_editor.iterrows():
            pos_in_base_editor = self.find_pos_in_base_editor(cds, row['pos'], row['targetSeq'], row['strand'])

            if self.killed_by_exon_filter(row, exon_boundaries, cds):
                continue
            elif self.killed_by_no_AS_mutation_ABE(row, pos_in_base_editor, cds):
                continue
            elif self.killed_by_no_AS_mutation_CBE(row, pos_in_base_editor, cds):
                continue
            else:
                unexpected.append(row)

        self.assertEqual(len(unexpected), 0, f"Unexplained guides: {len(unexpected)}{unexpected}")
        
