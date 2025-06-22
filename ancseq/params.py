import sys
import argparse
from ancseq.__init__ import __version__


class Params(object):

    def __init__(self, program_name):
        self.program_name = program_name

    def set_options(self):
        if self.program_name == 'ancseq':
            parser = self.ancseq_options()

        if len(sys.argv) == 1:
            args = parser.parse_args(['-h'])
        else:
            args = parser.parse_args()
        return args

    def ancseq_options(self):
        parser = argparse.ArgumentParser(description='ancseq 버전 {}'.format(__version__),
                                         formatter_class=argparse.RawTextHelpFormatter)
        parser.usage = 'ancseq -s <정렬된_FASTA> -m <모드> -o <출력_디렉토리> [-t <정수>]'

        # 옵션 설정
        parser.add_argument('-s',
                            '--seq',
                            action='store',
                            required=True,
                            type=str,
                            help='FASTA 형식의 서열 정렬 파일.',
                            metavar='')

        parser.add_argument('-m',
                            '--mode',
                            required=True,
                            type=str,
                            help='서열 유형. [DNA/AA/CODON]',
                            choices=['DNA','AA', 'CODON'],
                            metavar='')

        parser.add_argument('-o',
                            '--out',
                            action='store',
                            required=True,
                            type=str,
                            help='출력 디렉토리. 지정된 이름은 존재하지 않아야 합니다.',
                            metavar='')

        parser.add_argument('-t',
                            '--threads',
                            action='store',
                            default=4,
                            type=int,
                            help='스레드 수. [4]',
                            metavar='')
        
        parser.add_argument('-b',
                            '--bootstrap',
                            action='store',
                            default=1000,
                            type=int,
                            help='부트스트랩 반복 횟수. [1000]',
                            metavar='')
        
        parser.add_argument('--max-report',
                            action='store',
                            default=5,
                            type=int,
                            help='동일 위치에서 보고할 모호한 사이트의 최대 수. [5]',
                            metavar='')
        
        parser.add_argument('--min-prob',
                            action='store',
                            default=0.05,
                            type=float,
                            help='모호한 사이트로 보고될 최소 확률. [0.05]',
                            metavar='')
        
        parser.add_argument('--min-gap-prob',
                            action='store',
                            default=0.5,
                            type=float,
                            help='조상 상태를 갭으로 대체할 최소 확률. [0.5]',
                            metavar='')
        
        parser.add_argument('--fast',
                            action='store_true',
                            help='IQ-TREE에서 -fast 옵션 사용 [FALSE]')
        
        parser.add_argument('--model',
                            action='store',
                            default='MFP',
                            type=str,
                            help='IQ-TREE의 치환 모델 지정. 기본적으로 IQ-TREE는 ModelFinder를 사용하여 최적의 치환 모델을 검색합니다 [MFP]',
                            metavar='')

        parser.add_argument('--outgroup',
                            action='store',
                            default=None,
                            type=str,
                            help='IQ-TREE의 외부 그룹 지정. [None]',
                            metavar='')
        
        parser.add_argument('--stop-codon-prob',
                            action='store_true',
                            help='DNA 모드에서 코돈 확률 계산 중지 [FALSE]')
        
        parser.add_argument('--asr-only',
                            action='store_true',
                            help='계통수 구축을 건너뛰고 조상 상태만 재구축 [FALSE]')

        # 버전 설정
        parser.add_argument('-v',
                            '--version',
                            action='version',
                            version='%(prog)s {}'.format(__version__))

        return parser

    