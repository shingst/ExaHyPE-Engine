"""
.. module:: fileWriter
  :platform: Unix, Windows, Mac
  :synopsis: Provides routines to write vectors and matrices to file.
.. moduleauthor:: Angelika Schwarz <angelika.schwarz@tum.de>, Dominic Etienne Charrier <dominic.e.charrier@durham.ac.uk>

:synopsis: Provides routines to write vectors and matrices to file.
"""
def writeVectorToFile(vector, length, vectorname, filename):
    with open(filename,"a") as out:
        out.write("  const double exahype::dg::"+vectorname+"["+str(length)+"] = \n")
        out.write("  {\n")
        row = "    "
        for l in range(0,length):
            value = vector[l]
            string = "%s" %value
            row += string
            if l < length-1:
                row +=", \t"
        out.write(row)
        out.write("\n  };\n\n")
        out.close()
"""
   Writes a (n x m) matrix to a file using append mode.
"""
def writeMatrixToFile(matrix, n, m, matrixname, alignmentFlag, namespace, filename):
    
    if alignmentFlag :
        with open(filename,"a") as out:
            out.write("  const double "+namespace+matrixname+"["+str(n)+"]["+str(m)+"] "
                      " __attribute__((aligned(ALIGNMENT))) = \n")
            out.write("  {\n")
    else :
        with open(filename,"a") as out:
            out.write("  const double "+namespace+matrixname+"["+str(n)+"]["+str(m)+"]  = \n")
            out.write("  {\n")
        
    with open(filename,"a") as out:    
        for k in range(0,n):
            row = "    {"           
            for l in range(0,m):
                value = matrix[k][l]
                string = "%s" %value
                row += string
                if l < m-1:
                    row +=", \t"
               
            if k < n-1:
                out.write(row+"},\n")
            else:
                out.write(row+"}\n")
            
        out.write("  };\n\n")
        out.close()

"""
   Writes the entries of n vector to a file "out" 
"""
def writeVectorLookupTableInitToFile(out, indent, vector, n, vectorname):
    text = ""
    for l in range(0,n):
        value = vector[l]
        line = "%s%s[%d] = %s;\n" % (indent,vectorname,l,value)
        text += line
    out.write(text)

"""
   Writes the entries of (n x m) matrix to a file "out" 
"""
def writeMatrixLookupTableInitToFile(out, indent, matrix, n, m, matrixname):
    text = ""
    for k in range(0,n):
        for l in range(0,m):
            value = matrix[k][l]
            line = "%s%s[%d][%d] = %s;\n" % (indent,matrixname,k,l,value)
            text += line
    out.write(text)
