from __future__ import print_function
import re,os,time
import argparse
import shutil
import ipdb
from six import string_types as basestring

parser = argparse.ArgumentParser()
parser.add_argument("--reconvert",default=False,action='store_true')
#parser.add_argument("--arxiv",default=True,action='store_true')
parser.add_argument("--aanda",default=False,action='store_true')
parser.add_argument("--texit",default=False,action='store_true')
args = parser.parse_args()
print("ARGS: ",args)

ppath='/Users/adam/work/apex_cmz_h2co/tex/'
paper_name = 'apex_cmz_h2co'
file = open(os.path.join(ppath,paper_name+'.tex'),'r')

isarxiv = not args.aanda
outdir = "aanda/" if not isarxiv else "arxiv/"
outtype = "aanda" if not isarxiv else "arxiv"
print("Converting package to {0} format".format(outtype))

if not os.path.exists(ppath+outdir):
    os.mkdir(ppath+outdir)

shutil.copy('aamacros.tex', outdir)
shutil.copy('aa.cls', outdir)

outfn = ppath+outtype+'form_temp.tex'
print("Creating temporary file: {0}".format(outfn))
outf = open(outfn,'w')

inputre = re.compile('input{(.*)}')
includere = re.compile('include{(.*)}')
bibre = re.compile('bibliography{(.*)}')
aandare = re.compile("documentclass{aa}")
#                           \documentclass{aa}
beginre = re.compile('\\begin{document}')
endre = re.compile('\\end{document}')
prefacere = re.compile('\\input{preface.*}')
solobibre = re.compile("\\input{(solobib)}")

def strip_input(list_of_lines):
    # strip out preface, solobib, end{doc}
    #return "".join(list_of_lines[1:-2])
    return "".join(
                   [line
                    for line in list_of_lines
                    if not prefacere.search(line)
                    #and not solobibre.search(line)
                    and not endre.search(line)
                    and not beginre.search(line)
                   ]
    )

def dobib(bib, outf):
    #bn = bib.groups()[0] + '.bbl'
    bn = paper_name+".bbl"
    print("Doing bibliography " + bn)
    with open(ppath+bn,'r') as f:
        print(strip_input(f.readlines()), end='', file=outf)

spacecount=0
for ii,line in enumerate(file.readlines()):
    if line[0] == "%":
        continue
    elif line[0] == " ":
        spacecount += 1
    input = inputre.search(line)
    include = includere.search(line)
    bib = bibre.search(line)
    solobib = solobibre.search(line)
    if solobib is not None:
        fn = solobib.groups()[0] + ".tex"
        print(ii, "Doing solobib " + fn)
        with open(os.path.join(ppath,fn),'r') as f:
            solobib = f.readlines()
            for ln in solobib:
                if ln[0].strip() == "%":
                    continue
                bib = bibre.search(ln)
                if bib is not None:
                    print(ii,"Bib matched: ",bib.groups())
                    dobib(bib,outf)
    elif input is not None:
        fn = os.path.splitext(input.groups()[0])[0] + ".tex"
        if fn.count('.') > 1:
            ipdb.set_trace()
        print(ii, "Doing input " + fn)
        with open(os.path.join(ppath,fn),'r') as f:
            if 'preface' in line:
                print(f.read(), end='', file=outf)
            else:
                print(strip_input(f.readlines()), end='', file=outf)
    elif include is not None:
        fn = os.path.splitext(include.groups()[0])[0] + ".tex"
        print(ii, "Doing include " + fn)
        f = open(os.path.join(ppath,fn),'r')
        if 'preface' in line:
            print(f.read(), end='', file=outf)
        else:
            print(strip_input(f.readlines()), end='', file=outf)
        f.close()
    elif bib is not None:
        print(ii, "Doing bib (no solo) " + fn)
        dobib(bib,outf)
    else:
        print(line, end="", file=outf, sep="")

print("Spacecount = {0}".format(spacecount))

outf.close()
file.close()

file = open(outfn,'r')
outfn = ppath+outtype+'form.tex'
outf = open(outfn,'w')

default_suffix='.png'
figre = re.compile('Figure{{?(.*?)}?}')
#figre = re.compile('Figure\n?{{?(.*?)}?}?',flags=re.MULTILINE)
figre1line = re.compile('^\\\Figure',)
fig2re = re.compile('R?o?t?FigureTwoA?A?\n?{{?(.*?)}?}\n?\s*{(.*?)}?}',flags=re.MULTILINE)
fig2re1line = re.compile('^\\\R?o?t?FigureTwoA?A?\n?({{?(.*?)}?})?', flags=re.MULTILINE)
fig3re = re.compile('R?o?t?FigureThreeA?A?\n?{{?(.*?)}?}\n?\s*{(.*?)}?}',flags=re.MULTILINE)
fig3re1line = re.compile('^\\\R?o?t?FigureThreeA?A?\n?({{?(.*?)}?})?', flags=re.MULTILINE)
fig4re = re.compile('FigureFourP?D?F?\n?\s*{{?(.*?)}?}\n\s*?{{?(.*?)}}?\n\s*?{{?(.*?)}}?\n\s*?{{?(.*?)}}?',flags=re.MULTILINE)
fig4re1line = re.compile('FigureFourP?D?F?({{?(.*?)}?})?',flags=re.MULTILINE)
plotonere = re.compile('plotone{{?(.*?)}')
includegre = re.compile('includegraphics\[.*\]{{?(.*?)}?}')
fig_suffixes = "png|eps|pdf"
lonelygraphics = re.compile('^\s*{{?(.*?(%s)?)}?}' % fig_suffixes)

out_suffix = "pdf" if not isarxiv else "pdf"
out_suffix = 'pdf'

count = 1

for ii,line in enumerate(file.readlines()):
    aanda = aandare.search(line)
    if line[0] == "%":
        continue
    #elif aanda is not None:
    #    print "AandA -> aastex"
    #    print >>outf,"\\documentclass[12pt,preprint]{aastex}",
    #    f.close()
    #    continue
    fig1b = figre1line.search(line)
    fig1 = figre.search(line)
    fig2 = fig2re.search(line)
    fig2b = fig2re1line.search(line)
    fig3 = fig3re.search(line)
    fig3b = fig3re1line.search(line)
    fig4 = fig4re.search(line)
    fig4b = fig4re1line.search(line)
    pone = plotonere.search(line)
    lonely = lonelygraphics.search(line)
    igr = includegre.search(line)

    if igr is not None:
        fign = igr.groups()[0]
        prevline = ''
    elif fig1 is not None:
        fign = fig1.groups()[0]
        prevline = ''
    elif fig2 is not None:
        fign = fig2.groups()[0:2]
        prevline = ''
    elif fig3 is not None:
        fign = fig3.groups()[0:3]
        prevline = ''
    elif fig4 is not None:
        fign = fig4.groups()[0:4]
        prevline = ''
    elif fig4b is not None:
        fign = [n for n in fig4b.groups() if n is not None]
        nfigs = len(fign)
        prevline = 'fig4'
    elif fig2b is not None:
        fign = fig2b.groups()[1]
        if fign in (None,''):
            fign = None
            nfigs = 0
        else:
            nfigs = 1
        prevline = 'fig2'
    elif fig3b is not None:
        fign = fig3b.groups()[1]
        if fign in (None,''):
            fign = None
            nfigs = 0
        else:
            nfigs = 1
        prevline = 'fig3'
    elif pone is not None:
        fign = pone.groups()[0]
        prevline = ''
    elif fig1b is not None:
        if 'Two' in line:
            print("Two in line: ",line)
            ipdb.set_trace()
        #print "Found solo figure"
        fign = None#fig1b.groups()[0]
        nfigs = 0
        prevline = 'fig1'
    elif lonely is not None and prevline in ('fig4','fig2','fig3','fig1'):
        # DEBUG print "in lonely: ",lonely.groups()," nfigs=",nfigs," prevline=",prevline
        fign = lonely.groups()[0]
        nfigs += 1
        if nfigs >= int(prevline[-1]):
            prevline = ''
    elif lonely is not None:
        pass
        #print "found a lonely non-figure: ",lonely.group(), lonely.groups()
    else:
        fign=None
        #if 'prevline' in locals() and prevline != '':
        #    print "{1}: Found prevline={0}, setting to ''".format(prevline, line.strip())
        prevline = ''

    if line.strip() == "\FigureTwoAA":
        assert prevline == 'fig2'

    input = inputre.search(line)

    if fign is not None:
        #DEBUG print "Found fign: %s" % fign
        if fign in ("#1","#2","#3","#4"):
            print(line, end="", file=outf)
            continue
        if len(fign) == 1 or isinstance(fign,basestring):
            figlist = [fign]
        else:
            figlist = fign
        print("Figlist: ",figlist)
        outline = line
        for fign in figlist:
            fignroot = fign
            if fign[-4] != '.': # if no suffix, add one
                fign += default_suffix
            if fign[-3:] == "png":
                #if args.reconvert:
                #    os.system("pngtoeps %s" % fign)
                fn = ppath+fign.replace("png",out_suffix)
            elif fign[-3:] == "svg":
                #if args.reconvert:
                #    os.system("svg2eps %s" % fign)
                fn = ppath+fign.replace("svg",out_suffix)
            elif fign[-3:] == "pdf":
                #if args.reconvert:
                #    os.system("pdf2ps %s %s" % (fign,fign.replace("pdf","ps")))
                #    os.system("mv %s %s" % (fign.replace('pdf','ps'),fign.replace('pdf','eps')))
                fn = ppath+fign.replace("pdf",out_suffix)
            elif fign[-3:] != out_suffix:
                nroot = os.path.splitext(fign)[0]
                fn = os.path.join(ppath,nroot+"."+out_suffix)
            elif fign[-3:] == "tex":
                raise TypeError("Figure is a text file? %s" % fign)
            else:
                fn = ppath+fign
            if os.system('ls %s' % fn) != 0:
                ipdb.set_trace()
            print("Converting figure " + fn + " to f%i.%s" % (count,out_suffix))
            outfig = 'f%i.%s' % (count,out_suffix)
            outpath = os.path.join(ppath, outdir, outfig)
            if isarxiv:
                rslt = os.system('gs -dSAFER -dBATCH -dNOPAUSE -dAutoRotatePages=/None -dPDFSETTINGS=/screen -sDEVICE=pdfwrite -sOutputFile={1} {0}'.format(fn, outpath))
                if rslt != 0:
                    ipdb.set_trace()
            else:
                if os.system('cp %s %s' % (fn, outpath)) != 0:
                    ipdb.set_trace()
            if isarxiv:
                outline = outline.replace(fignroot,os.path.splitext(outfig)[0])
            else:
                #print(outline.replace(fignroot,outfig), " vs ", outline.replace(fignroot,os.path.splitext(outfig)[0]))
                outline = outline.replace(fignroot,outfig)
            count += 1
        print(outline, end="", file=outf)
    elif input is not None:
        fn = os.path.splitext(input.groups()[0])[0] + ".tex"
        print("Doing input " + fn)
        f = open(os.path.join(ppath,fn),'r')
        if 'preface' in line:
            print(f.read(), end="", file=outf)
        else:
            print(strip_input(f.readlines()), end="", file=outf)
        f.close()
    elif 'figures/' in line:
        raise ValueError("Found a line with figures/ in it, but it wasn't converted.")
    else:
        print(line, end="", file=outf)

outf.close()
file.close()

os.chdir(ppath)
os.system('cp %s %s/ms.tex' % (outfn,outdir))
if os.path.exists(outdir+'/Makefile'):
    os.chdir(outdir)
    os.system('make')
    os.system('rm ms.dvi')
    os.system('rm ms.ps')
    os.system('rm ms.aux')
    os.system('rm ms.log')
    os.chdir(ppath)
    os.system('mv %s/ms.pdf %s_draft%s.pdf' % (outdir,paper_name,time.strftime("%m%d",time.localtime())))
if args.texit:
    os.chdir(outdir)
    os.system('pdflatex ms.tex')
    os.system('pdflatex ms.tex')
    for junk in 'dvi','ps','aux','log':
        try:
            os.system('rm ms.'+junk)
        except IOError:
            pass
    os.chdir(ppath)
    if isarxiv:
        os.system('mv %s/ms.pdf %s_draft%s_arxiv.pdf' % (outdir,paper_name,time.strftime("%m%d",time.localtime())))
    else:
        os.system('mv %s/ms.pdf %s_draft%s_aanda.pdf' % (outdir,paper_name,time.strftime("%m%d",time.localtime())))

os.system('tar --exclude Makefile -czf Ginsburg_APEX_H2CO_CMZ_temmaps_%s_%s.tar.gz %s/ ' % (time.strftime("%m%d",time.localtime()),outtype,outdir))
