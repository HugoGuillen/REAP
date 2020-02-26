#Script takes as input several graphics and stores the summary to RNASeq_report.html

import sys, getopt, os, shutil, datetime
import argparse 


### parse options
parser = argparse.ArgumentParser(description="Creates a html report, storing all necessary files into html_summary. To visualize the report, keep the structure of the folder.")
parser.add_argument("-p", "--project", dest="project", help="name of the project")
parser.add_argument("-c", "--customer", dest="customer", help="group name of the customer", nargs="*")
parser.add_argument("-e", "--email", dest="email", help="email of the person who analyzed the data")
parser.add_argument("-n", "--names", dest="names", help="vector of all sample names", nargs="*")
parser.add_argument("-f", "--fastqc", dest="fastqc", help="vector of all sample names for the fastqc quality reports", nargs="*")
parser.add_argument("-u", "--unpaired", action="store_true", dest="unpaired", default=False, help="set if reads are unpaired (default: paired)") 
parser.add_argument("-q", "--quality", dest="quality", help="per base sequence quality graph(s) (fastqc): <Sample1 R1> <Sample1 R2> <Sample2 R1> <Sample2 R2> [optional]", nargs="*")
parser.add_argument("-g", "--gccont", dest="gccontent", help="GC-content graph (qualimap) [optional]")
parser.add_argument("-a", "--adapter", dest="adapter", help="adapter content graph(s) (fastqc): <Sample1 R1> <Sample1 R2> <Sample2 R1> <Sample2 R2> [optional]", nargs="*" )
parser.add_argument("-s", "--stat", dest="statistics", help="alignment statistics text files (qualimap, rnaseq_qc_results.txt) [optional]", nargs="*")
parser.add_argument("-d", "--dist", dest="distribution", help="read distribution statistics text files (Rseqc, readDistributionStats.txt) [optional]", nargs="*")
parser.add_argument("-b", "--geneb", dest="genebody", help="gene body coverage graph (Rseqc) [optional]", nargs="*")
parser.add_argument("-i", "--insert", dest="insertsize", help="insert size histogram graph (qualimap) [optional]")
parser.add_argument("-o", "--counts", dest="countspersample", help="counts per sample graph (qualimap) [optional]")
parser.add_argument("-t", "--sat", dest="saturation", help="saturation report graph(s) (Rseqc) [optional]", nargs="*")
parser.add_argument("-y", "--type", dest="genetype", help="counts per gene type graph(s) (qualimap) [optional]", nargs="*")
parser.add_argument("-z", "--pca", dest="pca", help="PCA plot (DESeq2) [optional]")
parser.add_argument("-m", "--ma", dest="ma", help="MA plot (DESeq2) [optional]")
parser.add_argument("-x", "--additional", dest="additional", help="Additional graph(s) you would like to append [optional]", nargs="*")


options = parser.parse_args()



### ------------------------ Functions -------------------------- ###


### create a labelled graph: the labelling is adjusted depending if the reads are unpaired or paired 
def labelled_figure(options_argument, options_names):
	HTML_ELEMENT = ""
	if options.unpaired == False:
		l = len(options_names)
                for i in range(0,l):
                        item1 = options_argument[(i+1)*2-2]
                        item2 = options_argument[(i+1)*2-1]
                        f1 = options_names[i] + "_R1_" + item1.split("/")[-1] #copy only filename without relative path
                        f2 = options_names[i] + "_R2_" + item2.split("/")[-1]
			fig_name1 = "Sample " + options_names[i] + " R1"
                        fig_name2 = "Sample " + options_names[i] + " R2"
			if item1 == "empty":
				s1 = (  '<figure class="block"><figure><figcaption>' + fig_name1 + '</figcaption><p>--- EMPTY PLOT ---</p></figure>\n\n')
			if item2 == "empty":
				s2 = (  '<figure class="block"><figure><figcaption>' + fig_name2 + '</figcaption><p>--- EMPTY PLOT ---</p></figure>\n\n')
			else:
				shutil.copy(item1, data_dir + f1)
        	                shutil.copy(item2, data_dir + f2)
                	        s1 = (  '<figure class="block"><figure><figcaption>' + fig_name1 + '</figcaption><img src="' + html_data + f1 + '" alt="' + fig_name1 + '"></figure>\n\n')
                        	s2 = (  '<figure><figcaption>' + fig_name2 + '</figcaption><img src="' + html_data + f2 + '" alt="' + fig_name2 + '"></figure></figure>\n\n')
                        HTML_ELEMENT += s1
                        HTML_ELEMENT += s2
	else:
		l = len(options_names)
                for i in range(0,l):
                        item = options_argument[i]
                        f = item.split("/")[-1] + options_names[i] #copy only filename without relative path
                        fig_name = "Sample " + options_names[i]
			if item == "empty":
				s = (  '<figure><figcaption>' + fig_name + '</figcaption>><p>--- EMPTY PLOT ---</p></figure>\n\n')
			else:
				shutil.copy(item, data_dir + f)
				s = (  '<figure><figcaption>' + fig_name + '</figcaption><img src="' + html_data + f + '" alt="' + fig_name + '"></figure>\n\n')
                        HTML_ELEMENT += s
	return HTML_ELEMENT

### create a simple figure
def simple_figure(options_argument):
	graph = options_argument
	f = graph.split("/")[-1]
        shutil.copy(graph, data_dir + f)
	if(is_img(options_argument)):
		s = '<img src="' + html_data + f + '" alt="' + f + '">\n\n'
       	else:
		s = '<embed src="' + html_data + f + '" alt="' + f + '">\n\n' 
	return s

### check format of input: return True if jpg or png, else False
def is_img(options_argument):
	last = options_argument[-3:]
	if(last == "jpg" or last == "png"):
		return True
	elif(last == "pdf"):
		return False
	else:
		print("\n\nWarning:\nImage format not jpg, png or pdf.\n\n")
		return False

### create simple labelled figure
def simp_fig(options_argument,label):
	graph = options_argument
        f = label + graph.split("/")[-1]
        shutil.copy(graph, data_dir + f)
	if(is_img(options_argument)):
		fig = '<figure><figcaption>' + label + '</figcaption><img src="' + html_data + f + '" alt="' + label + '"></figure>\n\n'
	else:
		fig = '<figure><figcaption>' + label + '</figcaption><embed src="' + html_data + f + '" alt="' + label + '"></figure>\n\n'
	return fig



def simple_labelled_figure(options_arguments): 
        i = 0
	figure = ""
        for graph in options_arguments:
                label = "Sample " + options.names[i]
                figure += simp_fig(graph, label)
                i += 1
	return figure



#################################################################################
#################################################################################




### store output and graphs in folder output
output_dir = "html_summary/"
html_data = "html_data/"
data_dir = output_dir+html_data

if not os.path.exists(data_dir):
	os.makedirs(data_dir)



### create the content of the html file
now = datetime.datetime.now()
HTML_INDEX = '''
	<html>\n
		<head>\n
			<meta charset="utf-8">
			<title>RNA-Seq report ''' + options.project + '''</title>\n
			<style>
				body {font-family:Tahoma, Geneva, sans-serif;}
				h1 {font-size:16px;}
				h2 {font-size:14px}
				h3, p, figcaption, table {font-size:12px;}
				img {max-width: 500px;}
				figure {margin: 2px; display:inline-block;}
				.block {min-width: 1010px; display:block;}
				table {display:inline-block; text-align:left; border-collapse:collapse; margin-right:5px;margin-left:0px;}
				td {border: 1px solid black; border-collapse: collapse; padding:5px;}
				th {border: 1px solid black; border-collapse: collapse; padding:5px;font-weight:bold;}
				caption {text-align:left; padding-top:10px; padding-bottom: 5px;}
				embed {width:500px;min-height:500px;}
			</style>\n
		</head>\n
		<body>\n
			<h1>RNA-seq report: Project ''' + options.project + '''</h1>\n
			
		<p>
		Group: ''' + " ".join(options.customer) + '''</br>
		Report created: '''+ now.strftime("%d %B %Y") + '''</br>
		Analysis for: ''' + options.email +'''</p>\n\n'''


CLOSING_TAG = '</body></html>'




#################################################################################
#################################################################################




### 1. Read quality report
READ_QUAL = '\n<h2>1. Read quality report</h2>\n' 
SEQ_QUAL = "\n\n"
GC_CONT = "\n\n"
ADP_CONT = "\n\n"


### 1.1 Per base sequence quality
if options.quality:
	SEQ_QUAL = '''<h3>1.1 Per base sequence quality (fastqc)</h3>\n
			<p>Base quality (Phred scores) along the length of the reads in each FastQ file. The box plots are drawn as follows: red line=median; yellow box=range between upper and lower quartiles; whiskers=range between 10% and 90% quantiles. The blue line shows the mean quality.</p>\n\n'''
	SEQ_QUAL += labelled_figure(options.quality, options.fastqc)
		
### 1.2 GC content
if options.gccontent:
	GC_CONT = '''<h3>1.2 GC content</h3>\n
			<p>Distribution of GC content of the reads for each sample. We would expect a roughly normal distribution centred on the average GC content of the genome which varies between species.</p>\n\n'''
	GC_CONT += simple_figure(options.gccontent)
	
### 1.3 Adapter content (fastqc)
if options.adapter:
	ADP_CONT = '''<h3>1.3 Adapter content (fastqc)</h3>\n
			<p>Cumulative proportion of observed adapter sequences along the length of the reads in each FastQ file. This should be (close to) zero.</p>\n\n'''
	ADP_CONT += labelled_figure(options.adapter, options.fastqc)
	
READ_QUAL += SEQ_QUAL + GC_CONT + ADP_CONT


### 2. Alignment report
ALIGN_REP = '\n<h2>2. Alignment report</h2>\n'
ALIGN_STAT = "\n\n"
READ_DIST = "\n\n"
GENE_BOD = "\n\n"
INS_SIZE= "\n\n"

### 2.1. Alignment statistics
### => Percent of reads mapped, total reads, unmapped reads
if options.statistics:
	ALIGN_STAT = '''<h3>2.1. Alignment statistics</h3>\n
		<p>Summary statistics for each sample:<br />
		i) Total reads=Total number of reads from the sequencer (for paired-end data, each pair counts as 2 reads)<br />
		ii) Number of mapped reads=Number of reads that can be mapped to the reference genome (for paired-end data, left=read 1, right=read 2)<br />
		iii) Non-unique alignments=Reads that match multiple positions in the reference genome<br />
		iv) Aligned to genes=Number of reads overlapping with an annotated gene<br />
		v) Ambiguous alignments=Number of reads that could derive from several genes<br />
		vi) No feature assigned=Reads that do not overlap with an annotated gene (intronic or intergenic)<br />
		vii) Not aligned=Reads that could not be mapped to the reference genome.</p>\n\n'''

	def percent(number, total):
		perc = float(number)*100/total
		perc = str(round(perc,2))
		return  " (" + perc + "%) "
	numbers = []
	row1 = "<tr><th></th>"
	row2 = "<tr><td>Total reads</td>"
	if options.unpaired == False:
		row3 = "<tr><td>Number of mapped reads (left)</td>"
		row4 = "<tr><td>Number of mapped reads (right)</td>"
	else:
		row3 = "<tr><td>Number of mapped reads</td>"
	row5 = "<tr><td>Number of non-unique alignments</td>"
	row6 = "<tr><td>Aligned to genes</td>"
	row7 = "<tr><td>Ambiguous alignments</td>"
	row8 = "<tr><td>No feature assigned</td>"
	row9 = "<tr><td>Not aligned</td>"
	
	i = 0 
	for element in options.statistics:
		stats = open(element,"r")
		stats = stats.readlines()
		if element == "empty":
			row1 += "<th>Sample " + options.names[i] + "</th>"
			row2 += "<td></td>"
			row3 += "<td></td>"
                        if options.unpaired == False:
				row4 += "<td></td>"
                        row5 += "<td></td>"
                        row6 += "<td></td>"
                        row7 += "<td></td>"
                        row8 += "<td></td>"
                        row9 += "<td></td>"
		else:
			if options.unpaired == False:
				left_right = stats[13].strip().split("=")[1]
				left = int(left_right.split("/")[0].replace(',',''))
				right = int(left_right.split("/")[1].replace(',',''))
				non_unique = int(stats[17].strip().split("=")[1].replace(',',''))
	                        to_genes = int(stats[18].strip().split("=")[1].replace(',',''))
        	                ambiguous = int(stats[19].strip().split("=")[1].replace(',',''))
                	        no_feature = int(stats[20].strip().split("=")[1].replace(',',''))
                        	not_aligned = int(stats[21].strip().split("=")[1].replace(',',''))

			else:
				left_right = int(stats[13].strip().split("=")[1].replace(',',''))
				non_unique = int(stats[16].strip().split("=")[1].replace(',',''))
				to_genes = int(stats[17].strip().split("=")[1].replace(',',''))
				ambiguous = int(stats[18].strip().split("=")[1].replace(',',''))
				no_feature = int(stats[19].strip().split("=")[1].replace(',',''))
				not_aligned = int(stats[20].strip().split("=")[1].replace(',',''))
		
			if options.unpaired == False:
				total = left + right + not_aligned
			else:
				total = left_right + not_aligned

		
			row1 += "<th>Sample " + options.names[i] + "</th>"
			row2 += "<td>" + str(format(total, ',d')) + "</td>"
			if options.unpaired == False:
				row3 += "<td>" + str(format(left,',d')) + percent(left,total) + "</td>"
				row4 += "<td>" + str(format(right,',d')) + percent(right,total) + "</td>"
			else:
				row3 += "<td>" + str(format(left_right,',d')) + percent(left_right,total) + "</td>"
			row5 += "<td>" + str(format(non_unique,',d')) + percent(non_unique,total)  + "</td>"
			row6 += "<td>" + str(format(to_genes,',d'))  + percent(to_genes,total) + "</td>"
			row7 += "<td>" + str(format(ambiguous,',d')) + percent(ambiguous,total)  + "</td>"
			row8 += "<td>" + str(format(no_feature,',d'))  + percent(no_feature,total) + "</td>"
			row9 += "<td>" + str(format(not_aligned,',d'))  + percent(not_aligned,total) + "</td>"
		i += 1		
	row1 += "</tr>"
	row2 += "</tr>"
	row3 += "</tr>"
	if options.unpaired == False:
		row4 += "</tr>"
	row5 += "</tr>"
	row6 += "</tr>"
	row7 += "</tr>"
	row8 += "</tr>"
	row9 += "</tr>"
	if options.unpaired == False:
		table = "\n\n<table>\n" + row1 + row2 + row3 + row4 + row5 + row6 + row7 + row8 + row9 + "\n</table>\n\n"
	else:
		table = "\n\n<table>\n" + row1 + row2 + row3 + row5 + row6 + row7 + row8 + row9 + "\n</table>\n\n"
	ALIGN_STAT += table

### 2.2 Read distribution 
if options.distribution:
	READ_DIST = '''<h3>2.2 Read distribution</h3>\n
			<p>Overview of how the reads are distributed among different genomic features (TSS=Transcription start site; TES=Transcription end site). Total_bases=Total length of each category in the genome; Tag_count=Number of reads assigned to each category. A read may be split and assigned to multiple categories; Tags/kb=Number of reads per kilobase of genomic sequence.</p>\n\n'''
	j = 0
	for file in options.distribution:
		table = '\n\n<table>\n<caption>Sample ' + options.names[j] + '</caption>'
		if file == "empty":
			table += '\n<tr><th></th><th></th><th></th><th></th></tr><tr><td></td><td></td><td></td><td></td></tr></table>\n\n'
		else:
			read_dist = open(file, "r")
			read_dist = read_dist.readlines()

			row = read_dist[4].split()
			th = '<tr><th>'+row[0]+'</th>'+'<th>'+row[1]+'</th>'+'<th>'+row[2]+'</th>'+'<th>'+row[3]+'</th></tr>'
			table += th
	
			for i in range(5,len(read_dist)-1):
				row = read_dist[i].split()
				tr = '<tr><td>'+row[0]+'</td>'+'<td>'+format(int(row[1]),',d')+'</td>'+'<td>'+format(int(row[2]),',d')+'</td>'+'<td>'+row[3]+'</td></tr>'
				table = table + tr
		table += '\n</table>\n\n'
		READ_DIST += table
		j += 1
	

### 2.3 Gene body coverage
if options.genebody:
	GENE_BOD = '''<h3>2.3 Gene body coverage</h3>\n
			<p>Distribution of reads along the length of the genes (5' end on left, 3' on right)</p>\n\n'''
	GENE_BOD += simple_labelled_figure(options.genebody)

### 2.4 Insert size histogram
if options.insertsize:
	INS_SIZE = '''<h3>2.4 Insert size histogram</h3>\n
			<p>Histogram of inferred insert size for each sample.</p>\n\n'''
	INS_SIZE += simple_figure(options.insertsize)

ALIGN_REP += ALIGN_STAT + READ_DIST + GENE_BOD + INS_SIZE


### 3. Counts report_REP = '<h2>3. Counts report</h2>'
COUNTS_REP = '<h2>3. Counts report</h2>\n'
COUNT_PER_SAMP = "\n\n"
SAT_REP = "\n\n"
COUNTS_PER_GENE = "\n\n"
PCA_PLOT = "\n\n"
MA_PLOT = "\n\n"

### 3.1. Counts per sample
if options.countspersample:
	COUNT_PER_SAMP = '''<h3>3.1. Counts per sample</h3>\n
				<p>Box plots of the distribution of counts across all genes in each sample.</p>\n\n'''
	COUNT_PER_SAMP += simple_figure(options.countspersample)

### 3.2 Saturation report
if options.saturation:
	SAT_REP = '''<h3>3.2 Saturation report</h3>\n
			<p>The number of splice junctions detected using different subsets of the data from 5% to 100% of all reads. Red=known junction based on the provided genome annotation; green=novel junctions; blue=all. At sequencing depths sufficient to perform alternative splicing analyses, at least the red line should reach a plateau where adding more data does not much increase the number of detected junctions.</p>\n\n'''
	SAT_REP += simple_labelled_figure(options.saturation)
		
### 3.3 Counts per gene type
if options.genetype:
	COUNTS_PER_GENE = '''<h3>3.3 Counts per gene type</h3>\n
				<p>Boxplots of the distribution of counts in each sample across all genes within a particular category.</p>\n\n'''
	COUNTS_PER_GENE += simple_labelled_figure(options.genetype)	

### 3.4 PCA
if options.pca:
	PCA_PLOT = '''<h3>3.4 PCA</h3>\n
			<p>Plot of the first two axes from a principal component analysis based on the 500 genes with the most variable expression across all samples. Different colours indicate the different experimental groups.</p>\n\n'''
	PCA_PLOT += simple_figure(options.pca)


### 3.5 MA plot
if options.ma:
	MA_PLOT = '''<h3>3.5 MA plot</h3>\n
			<p>Average number of reads observed per gene (mean expression) versus log2 fold-change. In red are genes showing differential expression using a threshold of P-adjusted < 5%.</p>\n\n'''
	MA_PLOT += simple_figure(options.ma)


COUNTS_REP += COUNT_PER_SAMP + SAT_REP + COUNTS_PER_GENE + PCA_PLOT + MA_PLOT



### write html
HTML_INDEX += READ_QUAL + ALIGN_REP + COUNTS_REP + CLOSING_TAG
outputfile = open(output_dir+"RNASeq_report.html", "w+")


outputfile.write(HTML_INDEX)








