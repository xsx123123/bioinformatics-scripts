#!/Users/zhangjian/miniconda3/bin/python
###### STAR ######
# AUTHOR  : zhang jian
# date    : 2024.1.17
# version : 2.1v
# loading model
######################
import sys
import os
import getopt
import logging
from datetime import datetime
from Bio import SeqIO
import logomaker as lm
######################
# define log_and_print function
def log_and_print(message):
    log_message = f'{datetime.now()} - {message}'
    logging.info(log_message)
    print(message)

# define run_logo function
def run_logo(text):
    print('\033[1;34m[INFO]\033[0m',datetime.now().strftime("%Y-%m-%d %H:%M:%S"),' \033[1;34m {}\033[0m'.format(text))

# define getparameter function
def getparameter(argv):
    if len(argv) == 0:  # 如果没有提供任何参数
        print(" \033[0;34m#------------------------------------------------------------------------------#\033[0m")
        print("")
        print("   \033[7;36mUSE :logo.py -f <fasta> -b <basenum>  -p <plottitle> -h <height> -w <width> \033[0m")
        print("")
        print(" \033[0;34m#------------------------------------------------------------------------------#\033[0m")
        sys.exit(2)
    fasta = ''
    basenum = ''
    plottitle = ''
    height = 2
    width = 6
    try:
        opts, args = getopt.getopt(argv, 'Hf:b:p:h:w:v', ['help', 'fasta=', 'basenum=',"plottitle=","height=","width=", 'version'])
    except getopt.GetoptError:
        print("logo.py -f <fasta> -b <basenum> -p <plottitle> -h <height> -w <width>")
        sys.exit(2)
    for opt, arg in opts:
        if opt in ('-H', '--help'):
            print(" \033[0;36m###############################################################\033[0m")
            print("")
            print(" \033[7;36m USE :logo.py -f <fasta> -b <basenum> -h <height> -w <width> \033[0m")
            print("")
            print("      \033[1;31m -f <fasta>\033[0m       :  😂 fasta file dir (Absolute Directory) ")
            print("      \033[1;31m -b <basenum>\033[0m     :  🤗 base cutoff ")
            print("      \033[1;31m -p <plottitle>\033[0m   :  😑 plot title ")
            print("      \033[1;31m -h <height>\033[0m      :  😟 plot height (dafult : 2) ")
            print("      \033[1;31m -w <width>\033[0m       :  😳 plot width (dafult : 6) ")
            print("      \033[1;31m -H <help>\033[0m        :  😳 print help ")
            print("")
            print(" \033[0;36m###############################################################\033[0m")
            sys.exit()
        elif opt in ('-f', '--fasta'):
            fasta = arg
        elif opt in ('-b', '--basenum'):
            basenum = arg
        elif opt in ('-p', '--plottitle'):
            plottitle = arg
        elif opt in ('-h', '--height'):
            height = arg
        elif opt in ('-w', '--width'):
            width = arg
        elif opt in ('-v', '--version'):
            sys.exit()
    return fasta,basenum,plottitle,height,width

# define run condition
def LOGORun(fasta, basenum,plottitle):
    directory = os.path.dirname(fasta)
    filename_without_extension = plottitle
    output_file_path_png = os.path.join(directory, f"{filename_without_extension}_logo.png")
    output_file_path_pdf = os.path.join(directory, f"{filename_without_extension}_logo.pdf")
    run_logo(("fasta file save dir : " + fasta))
    run_logo(("Logo plot(pdf) save dir : " + output_file_path_pdf))
    run_logo(("Logo plot(png) save dir : " + output_file_path_png))
    run_logo(("Logo plot title : " + plottitle))
    run_logo(("Draw logo plot basecutoff : " + basenum))
    return output_file_path_png,output_file_path_pdf

# define draw_plot function
def draw_plot(fasta,basenum,plottitle,output_file_path_png,output_file_path_pdf,height,width):
    sequences = []
    with open(fasta, 'r') as file:
        for record in SeqIO.parse(file, 'fasta'):
            sequence = str(record.seq[:int(basenum)])
            sequences.append(sequence)
    counts_mat = lm.alignment_to_matrix(sequences)
    counts_mat = lm.transform_matrix(counts_mat, from_type='counts', to_type='probability')
    logo = lm.Logo(counts_mat,
                color_scheme='classic')
    # set axes labels
    logo.ax.set_title(plottitle, fontsize=16)
    logo.ax.set_xlabel('Position (nt)',fontsize=14)
    logo.ax.set_ylabel("probability(%)", labelpad=1,fontsize=14)
    # Save plot
    fig = logo.ax.get_figure()
    fig.set_figwidth(int(width))
    fig.set_figheight(int(height))
    fig.savefig(output_file_path_png, format='png', dpi=300, bbox_inches='tight')
    fig.savefig(output_file_path_pdf, format='pdf', dpi=300, bbox_inches='tight')

# main function
def main():
    # set loggong result
    fasta,basenum,plottitle,height,width = getparameter(sys.argv[1:])
    run_logo(("Draw logo plot Do !!!"))
    output_file_path_png,output_file_path_pdf = LOGORun(fasta, basenum,plottitle)
    draw_plot(fasta,basenum,plottitle,output_file_path_png,output_file_path_pdf,height,width)
    run_logo(("Draw logo plot Done !!!"))

# run main()
if __name__ == "__main__":
    main()
###### END ######