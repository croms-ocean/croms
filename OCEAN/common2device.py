import os
import re
import sys
import argparse
import shutil

#order of .h below matters
#        python common2device.py ocean2d.h scalars.h ocean3d.h  grid.h   \
#                               coupling.h private_scratch.h mixing.h    \
#                               forces.h work.h ncscrum.h averages.h      \
#                               lmd_kpp.h climat.h nbq.h sources.h \
#                               wkb_wwave.h boundary.h tides.h bbl.h


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='*')
    parser.add_argument('-o', '--outfile', default='copy_to_devices.h')
    args = parser.parse_args()
    files = args.infile
    file_out = args.outfile

    
   # regex = re.compile("(real|integer|character)?\s*((?!\bdimension\b)(\w*))(\s*(\()(\s*\w*\s*\,*){1,6}(\s*\)))",re.IGNORECASE)
    
  #  re_tab77 = re.compile(r"((?<=real)|(?<=integer)|(?<=logical)|(?<=character)).*(\w+)(?=((\s*\()(((\s*\w*\s*(-|\+|\*)?\s*)+\:?){1,2},?){1,6}(\s*\))))(?<!dimension|parameter)",re.IGNORECASE)
   
  #  re_tab77 = re.compile(r"(\w+)(?=((\s*\()(((\s*\w*\s*(-|\+|\*)?\s*)+\:?){1,2},?){1,6}(\s*\))))(?<!dimension|parameter)",re.IGNORECASE)

    re_tab77 = re.compile(r"((\w+?)(?=(\s*?\()((\s*?\w+?\s*?(-|\+|\*)?:?){0,2},?){0,6}(\s*?\))))(?<!dimension|parameter)",re.IGNORECASE)

    re_tab90 = re.compile(r"(?!.*parameter.*)(real|logical|integer|character)(.*dimension.*::\s*)(.*)\!?",re.IGNORECASE)

    re_comment = re.compile(r"(^\s*!)|(^C)",re.IGNORECASE)
    re_cpp = re.compile(r"(.*defined.*)|(^#)")

   # re_tab77 = re.compile(r"(?<=\breal\s)(\w+)")

    buffall = []
  #  firstline = '! '+file_in+'\n!$acc enter data if(compute_on_device) copyin('
    firstline = '!$acc enter data if(compute_on_device) copyin('
    lastline = '\n!$acc& )\n'
    foundit=0

    for file_in in files:
        file = open(file_in, 'r', encoding='UTF-8')
        buffline = ['\n!'+file_in]
        
        # Main loop
        for line in file:
            if re_comment.match(line):
                continue
            elif re_cpp.match(line):
                buffline.append((line.rstrip('\n')))
            else:
                listmot = ['!$acc&']
                i=1
                for match in re_tab77.finditer(line):
                    #print('all :',match.group(0))
                    if foundit == 0 and i==1:
                        foundit=1
                        listmot.append(' '+match.group(0)) 
                    else:
                        listmot.append(', '+match.group(0))
                for match in re_tab90.finditer(line):
                    #print('all :',match.group())
                    if foundit == 0 and i==1:
                        foundit=1
                        listmot.append(' '+match.group(3)) 
                    else:
                        listmot.append(', '+match.group(3))
                
                    i=i+1    
                if(len(listmot) > 1):
                    buffline.append(''.join(listmot))
        file.close()
        if foundit == 0:
            continue
        # Dirty cleaning
        i = 0
        for line in buffline:
            if len(line) > 0 and line[0] == '#':
                if re.match('#\s*if', line) \
                and re.match('#\s*endif', buffline[i+1]):
                    buffline[i] = ''
                    buffline[i+1] = ''
            i = i + 1

        buffall.extend(buffline)

# Output
    buffall.insert(0,firstline)
    buffall.append(lastline) 
    buffall = [i for i in buffall if i]
    with open((file_out), 'w') as file:
        file.write("\n".join(buffall))
