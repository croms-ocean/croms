import os
import re
import sys
import argparse
import tempfile
import shutil

def doloop2d_treat(result,lines,file_out_main):

#	print('lines = ',lines)
	(lines,cd1d) = re.subn(r"CD\(i,([^)]*)\)",r"CD1D(\1)",lines)
	(lines,dc1d) = re.subn(r"DC\(i,([^)]*)\)",r"DC1D(\1)",lines)
	(lines,cf1d) = re.subn(r"CF\(i,([^)]*)\)",r"CF1D(\1)",lines)
	(lines,ffc1d)= re.subn(r"FFC\(i,([^)]*)\)",r"FFC1D(\1)",lines)
	(lines,fc1d) = re.subn(r"FC\(i,([^)]*)\)",r"FC1D(\1)",lines)
	(lines,dz1d) = re.subn(r"dZ\(i,([^)]*)\)",r"dZ1D(\1)",lines)
	(lines,dr1d) = re.subn(r"dR\(i,([^)]*)\)",r"dR1D(\1)",lines)
	(lines,bc1d) = re.subn(r"BC\(i,([^)]*)\)",r"BC1D(\1)",lines)
	file_out_main.write("      DO j=%s,%s\n" % (result.group(3),result.group(4)))
	accprivate="!$acc loop private("
	comma=""
	if (cd1d>0):
		accprivate=accprivate+"CD1D"
		comma=","
	if (dc1d>0):
		accprivate=accprivate+comma+"DC1D"
		comma=","
	if (cf1d>0):
		accprivate=accprivate+comma+"CF1D"
		comma=","
	if (ffc1d>0):
		accprivate=accprivate+comma+"FFC1D"
		comma=","
	if (fc1d>0):
		accprivate=accprivate+comma+"FC1D"
		comma=","
	if (dz1d>0):
		accprivate=accprivate+comma+"dZ1D"
		comma=","
	if (dr1d>0):
		accprivate=accprivate+comma+"dR1D"
		comma=","
	if (bc1d>0):
		accprivate=accprivate+comma+"BC1D"
		comma=","		
	accprivate=accprivate+") vector\n"
	if (cd1d+dc1d+cf1d+ffc1d+fc1d+dz1d+dr1d+bc1d>0):
		file_out_main.write(accprivate)
	file_out_main.write("        DO i=%s,%s\n" % (result.group(1),result.group(2)))
	get_nextline=re.compile(r".*(\n|$)")
	doloop1d=r"^[ \t]*do[ \t]*i[ \t]*=[ \t]*%s[ \t]*,[ \t]*%s(.*)\n" % (result.group(1),result.group(2))
	enddoloop1d=r"^[ \t]*enddo(.*)\n"
	doloop1ddeclare=re.compile(doloop1d,re.IGNORECASE)
	enddoloop1ddeclare=re.compile(enddoloop1d,re.IGNORECASE)
	indoloop=False
	while lines:
		getloop1d = doloop1ddeclare.match(lines)
		if getloop1d:
			ibegin=getloop1d.end()
			indoloop=True
		else:
			if indoloop:
				getendloop1d = enddoloop1ddeclare.match(lines)
				if getendloop1d:
					ibegin=getendloop1d.end()
					indoloop=False
				else:
					nextline=get_nextline.match(lines)
					file_out_main.write("%s" % nextline.group(0))
					ibegin=nextline.end()
			else:
				nextline=get_nextline.match(lines)
				file_out_main.write("%s" % nextline.group(0))
				ibegin=nextline.end()
		lines=lines[ibegin:len(lines)]
	file_out_main.write("      ENDDO\n")
	file_out_main.write("      ENDDO\n")

def doextend_treat(result,lines,file_out_main):
	file_out_main.write("      DO %s=%s,%s\n" % (result.group(1),result.group(2),result.group(3)))
	lines = re.sub(r"FX\(([^)]*),([^)]*)\)",r"FX_3D(\1,\2,k)",lines)
	arrays=result.group(4).split(',')
	for array in arrays:
		replace_string_in=r"%s\(([^)]*),([^)]*)\)"  % (array)
		replace_string_out=r"%s_3D(\1,\2,k)" % (array)
		lines = re.sub(replace_string_in,replace_string_out,lines)
		#print("replace = ",replace_string_in,replace_string_out)
	file_out_main.write(lines)
	file_out_main.write("      ENDDO\n")

#	print('lines = ',lines)
	(lines,cd1d) = re.subn(r"CD\(i,([^)]*)\)",r"CD1D(\1)",lines)
if __name__=="__main__":
	parser=argparse.ArgumentParser()
	parser.add_argument('infile')
	parser.add_argument('outfile')
	args=parser.parse_args()
	file = open(args.infile,'r',encoding='iso-8859-1')
	file_out_main=open(args.outfile,'w',encoding='iso-8859-1')
	insubroutine=False
	listoflines=[]
	listofuses=[]
	myfile=''.join(file.readlines())
	doloop2ddeclare=re.compile(r"^[ \t]*DOLOOP2D\((.*),(.*),(.*),(.*)\).*\n")
	enddoloop2ddecalre=re.compile(r"[ \t]*ENDDOLOOP2D.*\n")
	doextenddeclare=re.compile(r"^[ \t]*DOEXTEND\(([^,]*),([^,])*,([^,])*,(.*)\).*\n")
	enddoextenddeclare=re.compile(r"[ \t]*ENDDOEXTEND.*\n")
	get_nextline=re.compile(r".*(\n|$)")
	while myfile:
		doloop2d = doloop2ddeclare.match(myfile)
		if doloop2d:
			ibegin=doloop2d.end()
			myfile=myfile[ibegin:len(myfile)]
			lines=""
			while myfile:
				enddoloop2d = enddoloop2ddecalre.match(myfile)
				if enddoloop2d:
					doloop2d_treat(doloop2d,lines,file_out_main)
					ibegin=enddoloop2d.end()
					myfile=myfile[ibegin:len(myfile)]
					break
				else:
					nextline=get_nextline.match(myfile)
					lines=lines+nextline.group(0)
					ibegin=nextline.end()
					myfile=myfile[ibegin:len(myfile)]
		else:
                    doextend=doextenddeclare.match(myfile)
                    if doextend:
                        ibegin=doextend.end()
                        myfile=myfile[ibegin:len(myfile)]
                        lines=""
                        while myfile:
                            enddoextend=enddoextenddeclare.match(myfile)
                            if enddoextend:
                                doextend_treat(doextend,lines,file_out_main)
                                ibegin=enddoextend.end()
                                myfile=myfile[ibegin:len(myfile)]
                                break
                            else:
                                nextline=get_nextline.match(myfile)
                                lines=lines+nextline.group(0)
                                ibegin=nextline.end()
                                myfile=myfile[ibegin:len(myfile)]
                    else:
                        nextline=get_nextline.match(myfile)
                        file_out_main.write("%s" % nextline.group(0))
                        ibegin=nextline.end()
                        myfile=myfile[ibegin:len(myfile)]

	file.close()
	file_out_main.close()
