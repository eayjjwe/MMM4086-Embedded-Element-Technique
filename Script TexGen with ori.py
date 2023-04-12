#eayjjwe@nottingham.ac.uk MMME4086 MEng Individual Project, "Embedded Element Technique: Enabling the Effective Meshing of 3D Composites"
#REQUIREMENTS: Mat0 to be matrix (iso)
#              Mat1 to be weave (ortho)
#              Matrix to be solid part
#              1x .inp and 1x .ori file to be provided for each of matrix and weave
#              Matrix to be exported from TexGen as ABAQUS Voxel mesh, with NO boundary conditions and yarns deleted

#Note: throughout comments/script, terms "Host" and "Matrix" are used interchangeably

#import numpy to allow for easy matrix operations in Material Correction section
import numpy as np
from numpy.linalg import inv

#INP file formatting -> f_1 = embedded; f_2 = matrix; f_3 = output; f_1_ori = embedded .ori file; f_2_ori = matrix .ori file; f_4 = intermediate .ori file;
#f_5 = intermediate .inp file
EE_in = "weave_only.inp"; #f_1
H_in = "matrix_only2.inp"; #f_2

INP_out = "max_node_test.inp"; #f_3

EE_ori = "weave_only.ori" #f_1_ori
H_ori = "matrix_only.ori" #f_2_ori

Int_ori = "IntermediateOri.ori" #f_4 - holds the new .ori file before it is written back to EE_ori or H_ori
Int_inp = "IntermediateINP.inp" #f_5 - holds the corrected EE nodal definitions before it is written back to INP_out

#Adding prefix of "Part_name." before node numbers in EE .ori file
with open(EE_ori, "r") as f_1_ori, open(Int_ori, "w") as f_4:
    for line in f_1_ori:
        if not line.startswith("*"):
            if line.startswith(","):
                f_4.write("EE_only." + line.strip(", ").lstrip())
                line = f_1_ori.readline()
            elif not line.startswith(","):
                f_4.write("EE_only." + line)
        else:
            f_4.write(line)

with open(EE_ori, "w") as f_1_ori:
    f_1_ori.write("")

with open(EE_ori, "a") as f_1_ori, open(Int_ori, "r") as f_4:
    for line in f_4:
        f_1_ori.write(line)

#Adding prefix of "Part_name." before node numbers in matrix .ori file
with open(H_ori, "r") as f_2_ori, open(Int_ori, "w") as f_4:
    for line in f_2_ori:
        if not line.startswith("*"):
            if line.startswith(","):
                f_4.write("Host_only." + line.strip(", ").lstrip())
                line = f_2_ori.readline()
            elif not line.startswith(","):
                f_4.write("Host_only." + line)
        else:
            f_4.write(line)

with open(H_ori, "w") as f_2_ori:
    f_2_ori.write("")

with open(H_ori, "a") as f_2_ori, open(Int_ori, "r") as f_4:
    for line in f_4:
        f_2_ori.write(line)

#Standard Heading
with open(INP_out,"w") as f_3:
    f_3.write("**EET INP file generated through Script TexGen WiP.py\n*Heading\n*Preprint, echo=NO, model=NO, history=NO, contact=NO\n")

#Correcting node positions for outermost nodes (so all nodes lie within EE)
    #Define maximum & minimum coords in 1,2,3 directions for the matrix
node="*Node"
end_node="*"

H_1=[]
H_2=[]
H_3=[]

with open(H_in,"r") as f_2:
    node_section_started=False
    for line in f_2:
        if node in line:
            node_section_started=True
            line = f_2.readline()
        elif node_section_started and end_node in line:
            break
        if node_section_started:
            H_1.append(float((line.strip().split(",")[1]).strip()))
            H_2.append(float((line.strip().split(",")[2]).strip()))
            H_3.append(float((line.strip().split(",")[-1]).strip()))

H_max1=max(H_1)
H_max2=max(H_2)
H_max3=max(H_3)

H_min1=min(H_1)
H_min2=min(H_2)
H_min3=min(H_3)

    #Define maximum & minimum coords in 1,2,3 directions for the matrix
node="*Node"
end_node="*"

EE_1=[]
EE_2=[]
EE_3=[]

with open(EE_in,"r") as f_1:
    node_section_started=False
    for line in f_1:
        if node in line:
            node_section_started=True
            line = f_1.readline()
        elif node_section_started and end_node in line:
            break
        if node_section_started:
            EE_1.append(float((line.strip().split(",")[1]).strip()))
            EE_2.append(float((line.strip().split(",")[2]).strip()))
            EE_3.append(float((line.strip().split(",")[-1]).strip()))

EE_max1=max(EE_1)
EE_max2=max(EE_2)
EE_max3=max(EE_3)

EE_min1=min(EE_1)
EE_min2=min(EE_2)
EE_min3=min(EE_3)

    #Redefining EE part coordinates which exceed the matrix coordinates
node="*Node"
end_node="*"

with open(EE_in, "r") as f_3, open(Int_inp, "w") as f_5:
    node_section_started=False
    for line in f_3:
        if node in line:
            node_section_started=True
            line=f_3.readline()
        elif node_section_started and end_node in line:
            break
        if node_section_started:
            if float((line.strip().split(",")[1]).strip()) > H_max1:
                f_5.write(str(int((line.strip().split(",")[0]).strip())) + ", " + str(H_max1) + ", " + str(float((line.strip().split(",")[2]).strip())) + ", " + str(float((line.strip().split(",")[-1]).strip())))
                f_5.write("\n")
            elif float((line.strip().split(",")[2]).strip()) > H_max2:
                f_5.write(str(int((line.strip().split(",")[0]).strip())) + ", " + str(float((line.strip().split(",")[1]).strip())) + ", " + str(H_max2) + ", " + str(float((line.strip().split(",")[-1]).strip())))
                f_5.write("\n")
            elif float((line.strip().split(",")[-1]).strip()) > H_max3:
                f_5.write(str(int((line.strip().split(",")[0]).strip())) + ", " + str(float((line.strip().split(",")[1]).strip())) + ", " + str(float((line.strip().split(",")[2]).strip())) + ", " + str(H_max3))
                f_5.write("\n")
            elif float((line.strip().split(",")[1]).strip()) < H_min1:
                f_5.write(str(int((line.strip().split(",")[0]).strip())) + ", " + str(H_min1) + ", " + str(float((line.strip().split(",")[2]).strip())) + ", " + str(float((line.strip().split(",")[-1]).strip())))
                f_5.write("\n")
            elif float((line.strip().split(",")[2]).strip()) < H_min2:
                f_5.write(str(int((line.strip().split(",")[0]).strip())) + ", " + str(float((line.strip().split(",")[1]).strip())) + ", " + str(H_min2) + ", " + str(float((line.strip().split(",")[-1]).strip())))
                f_5.write("\n")
            elif float((line.strip().split(",")[-1]).strip()) < H_min3:
                f_5.write(str(int((line.strip().split(",")[0]).strip())) + ", " + str(float((line.strip().split(",")[1]).strip())) + ", " + str(float((line.strip().split(",")[2]).strip())) + ", " + str(H_min3))
                f_5.write("\n")
            else:
                f_5.write(line)

#Creating node sets (NSets) of the extreme faces of the matrix body
H_nset_max1=[]
H_nset_min1=[]
H_nset_max2=[]
H_nset_min2=[]
H_nset_max3=[]
H_nset_min3=[]

node="*Node"
end_node="*"

with open(H_in, "r") as f_2:
    node_section_started=False
    for line in f_2:
        if node in line:
            node_section_started=True
            line=f_2.readline()
        elif node_section_started and end_node in line:
            break
        if node_section_started:
            #max 1
            if float((line.strip().split(",")[1]).strip()) == H_max1:
                H_nset_max1.append(int((line.strip().split(",")[0]).strip()))
            #min 1
            if float((line.strip().split(",")[1]).strip()) == H_min1:
                H_nset_min1.append(int((line.strip().split(",")[0]).strip()))
            #max 2
            if float((line.strip().split(",")[2]).strip()) == H_max2:
                H_nset_max2.append(int((line.strip().split(",")[0]).strip()))
            #min 2
            if float((line.strip().split(",")[2]).strip()) == H_min2:
                H_nset_min2.append(int((line.strip().split(",")[0]).strip()))            
            #max 3
            if float((line.strip().split(",")[-1]).strip()) == H_max3:
                H_nset_max3.append(int((line.strip().split(",")[0]).strip()))
            #min 3
            if float((line.strip().split(",")[-1]).strip()) == H_min3:
                H_nset_min3.append(int((line.strip().split(",")[0]).strip()))

    #formatting the lists to remove the [] at the start/end
H_nset_max1=str(H_nset_max1)
H_nset_max1=H_nset_max1[1:-1]

H_nset_min1=str(H_nset_min1)
H_nset_min1=H_nset_min1[1:-1]

H_nset_max2=str(H_nset_max2)
H_nset_max2=H_nset_max2[1:-1]

H_nset_min2=str(H_nset_min2)
H_nset_min2=H_nset_min2[1:-1]

H_nset_max3=str(H_nset_max3)
H_nset_max3=H_nset_max3[1:-1]

H_nset_min3=str(H_nset_min3)
H_nset_min3=H_nset_min3[1:-1]

#Part 1 (embedded element) from INP1
with open(INP_out,"a") as f_3:
    f_3.write("**\n**PARTS\n**\n*Part, name=EE_only\n")

    #Corrected node data from f_5
with open(INP_out,"a") as f_3, open(Int_inp,"r") as f_5:
    f_3.write("*Node\n")
    for line in f_5:
        f_3.write(line)    

    #Copying *Element to ***Materials*** from INP1, excluding orientation data (between orientation and element set)
Part = "*Element"
End_Part = "*** MATERIALS ***"

Dist_start = "*Distribution Table,"
Dist_end = "*Distribution,"

with open(EE_in, "r") as f_1, open(INP_out, "a") as f_3:
    node_section_started = False
    exclude_section_started = False
    for line in f_1:
        if Part in line:
            node_section_started = True
        elif node_section_started and line.strip() == End_Part:
            break
        elif Dist_start in line:
           exclude_section_started = True
        elif Dist_end in line:
           exclude_section_started = False
           line = f_1.readline()
        if node_section_started and not exclude_section_started:
            f_3.write(line)

    #copying solid section parts
Solid_Section = "*Solid Section"
Section_controls = "*Section Controls"

with open(EE_in, "r") as f_1, open(INP_out, "a") as f_3:
    solid_section_started = False
    for line in f_1:
        if Solid_Section in line:
            solid_section_started = True
        elif solid_section_started and Section_controls in line:
            break
        if solid_section_started:
            f_3.write(line)

with open(INP_out,"a") as f_3:
    f_3.write("*End Part\n")

#Part 2 (Host) from INP 2
with open(INP_out,"a") as f_3:
    f_3.write("**\n*Part, name=Host_only\n")

    #copying node to ***Materials*** from INP1
Part = "*Node"
End_Part = "*** MATERIALS ***"

Dist_start = "*Distribution Table,"
Dist_end = "*Distribution,"

with open(H_in, "r") as f_2, open(INP_out, "a") as f_3:
    node_section_started = False
    exclude_section_started = False
    for line in f_2:
        if Part in line:
            node_section_started = True
        elif node_section_started and line.strip() == End_Part:
            break
        elif Dist_start in line:
           exclude_section_started = True
        elif Dist_end in line:
           exclude_section_started = False
           line = f_2.readline()
        if node_section_started and not exclude_section_started:
            f_3.write(line)            

    #copying solid section parts - Mat0
SSolid_Section = "*Solid Section"

SS_found = True

with open(H_in, "r") as f_2, open(INP_out, "a") as f_3:
    for line in f_2:
        if Solid_Section in line:
            SS_found = False
            f_3.write(line)
            line = f_2.readline()
            f_3.write(line)
            break
        if SS_found:
            f_3.write("*Solid Section, Elset=Matrix, Material=Mat0, Orientation=TexGenOrientations, controls=HourglassEnhanced\n1.0,\n")
        break 
with open(INP_out,"a") as f_3:
    f_3.write("*End Part\n")

#Assembly Definition
with open(INP_out,"a") as f_3:
    f_3.write("**\n**ASSEMBLY\n**\n*Assembly, name=composite\n*Instance, name=EE_only, part=EE_only\n*End instance\n*Instance, name=Host_only, part=Host_only\n*End instance\n")

    #Define max nodes
Element_search = "*Element"

with open(EE_in, "r") as f_1:
    previous_line = None
    for line in f_1:
        if Element_search in line:
            if previous_line is not None:
                max_Node_EE = previous_line.split(",")[0].lstrip()
            break
        previous_line = line

with open(H_in, "r") as f_2:
    previous_line = None
    for line in f_2:
        if Element_search in line:
            if previous_line is not None:
                max_Node_H = previous_line.split(",")[0].lstrip()
            break
        previous_line = line

    #Define Nset & Elset for EE
    #Nset
with open(INP_out,"a") as f_3:
    f_3.write("*Nset, nset=EE_nset, instance=EE_only, generate\n1, " + max_Node_EE + ", 1\n")

    #Elset
Elset_EE = "*ElSet, ElSet=All"

with open(EE_in, "r") as f_1, open(INP_out, "a") as f_3:
    for line in f_1:
        if Elset_EE in line:
            f_3.write("*ElSet, ElSet=All_EE, instance=EE_only, generate\n")
            line = f_1.readline()
            f_3.write(line)
            break
    
    #Define Nset & Elset for Host
    #Nset
with open(INP_out, "a") as f_3:
    f_3.write("*Nset, nset=H_nset, instance=Host_only, generate\n1, " + max_Node_H + ", 1\n")

    #Elset
Elset_H1 = "*ElSet"
Elset_H2 = "Generate"

with open(H_in, "r") as f_2, open(INP_out, "a") as f_3:
    for line in f_2:
        if Elset_H1 and Elset_H2 in line:
            f_3.write("*ElSet, ElSet=Matrix, instance=Host_only, generate\n")
            line = f_2.readline()
            f_3.write(line)
            break

    #NSet generation for extreme matrix faces (max/min for 1,2,3 directions)
line_length = 40     #up to 16 entries per line -> min length (nodes 1-16) is ~40-50 depending if spaces are included

    #Max1
with open(INP_out, "a") as f_3:
    f_3.write("*Nset, nset=H_nset_max1, instance=Host_only\n")
    # loop through the input string
    while len(H_nset_max1) > 0:
        if len(H_nset_max1) <= line_length:
            f_3.write(H_nset_max1.strip())
            break
        # find the position of the closest comma before the line length
        comma_pos = H_nset_max1[:line_length+1].rfind(",")
        if comma_pos == -1:
            comma_pos = line_length
        # write the line up to the comma position
        f_3.write(H_nset_max1[:comma_pos+1].strip()[:-1] + "\n")
        H_nset_max1 = H_nset_max1[comma_pos+1:]

    #Min1
with open(INP_out, "a") as f_3:
    f_3.write("\n*Nset, nset=H_nset_min1, instance=Host_only\n")
    # loop through the input string
    while len(H_nset_min1) > 0:
        if len(H_nset_min1) <= line_length:
            f_3.write(H_nset_min1.strip())
            break
        # find the position of the closest comma before the line length
        comma_pos = H_nset_min1[:line_length+1].rfind(",")
        if comma_pos == -1:
            comma_pos = line_length
        # write the line up to the comma position
        f_3.write(H_nset_min1[:comma_pos+1].strip()[:-1] + "\n")
        H_nset_min1 = H_nset_min1[comma_pos+1:]

    #Max2
with open(INP_out, "a") as f_3:
    f_3.write("\n*Nset, nset=H_nset_max2, instance=Host_only\n")
    # loop through the input string
    while len(H_nset_max2) > 0:
        if len(H_nset_max2) <= line_length:
            f_3.write(H_nset_max2.strip())
            break
        # find the position of the closest comma before the line length
        comma_pos = H_nset_max2[:line_length+1].rfind(",")
        if comma_pos == -1:
            comma_pos = line_length
        # write the line up to the comma position
        f_3.write(H_nset_max2[:comma_pos+1].strip()[:-1] + "\n")
        H_nset_max2 = H_nset_max2[comma_pos+1:]

    #Min2
with open(INP_out, "a") as f_3:
    f_3.write("\n*Nset, nset=H_nset_min2, instance=Host_only\n")
    # loop through the input string
    while len(H_nset_min2) > 0:
        if len(H_nset_min2) <= line_length:
            f_3.write(H_nset_max1.strip())
            break
        # find the position of the closest comma before the line length
        comma_pos = H_nset_min2[:line_length+1].rfind(",")
        if comma_pos == -1:
            comma_pos = line_length
        # write the line up to the comma position
        f_3.write(H_nset_min2[:comma_pos+1].strip()[:-1] + "\n")
        H_nset_min2 = H_nset_min2[comma_pos+1:]

    #Max3
with open(INP_out, "a") as f_3:
    f_3.write("\n*Nset, nset=H_nset_max3, instance=Host_only\n")
    # loop through the input string
    while len(H_nset_max3) > 0:
        if len(H_nset_max3) <= line_length:
            f_3.write(H_nset_max3.strip())
            break
        # find the position of the closest comma before the line length
        comma_pos = H_nset_max3[:line_length+1].rfind(",")
        if comma_pos == -1:
            comma_pos = line_length
        # write the line up to the comma position
        f_3.write(H_nset_max3[:comma_pos+1].strip()[:-1] + "\n")
        H_nset_max3 = H_nset_max3[comma_pos+1:]

    #Min3
with open(INP_out, "a") as f_3:
    f_3.write("\n*Nset, nset=H_nset_min3, instance=Host_only\n")
    # loop through the input string
    while len(H_nset_min3) > 0:
        if len(H_nset_min3) <= line_length:
            f_3.write(H_nset_min3.strip())
            break
        # find the position of the closest comma before the line length
        comma_pos = H_nset_min3[:line_length+1].rfind(",")
        if comma_pos == -1:
            comma_pos = line_length
        # write the line up to the comma position
        f_3.write(H_nset_min3[:comma_pos+1].strip()[:-1] + "\n")
        H_nset_min3 = H_nset_min3[comma_pos+1:]

    #EE Constraint
with open(INP_out,"a") as f_3:
    f_3.write("\n**\n**CONSTRAINT - EET\n*Embedded Element, host elset=Matrix\nAll_EE\n")

    #End assembly
with open(INP_out,"a") as f_3:
    f_3.write("*End Assembly\n**\n")

#*Distribution Section
Part = "*Node"
End_Part = "*** MATERIALS ***"

Dist_start = "*Distribution Table,"
Dist_end = "*Distribution,"

with open(EE_in, "r") as f_1, open(INP_out, "a") as f_3:
    DistEE_section_started = False
    for line in f_1:
        if Dist_start in line:
            DistEE_section_started = True
        elif DistEE_section_started and Dist_end in line:
            f_3.write(line)
            break
        if DistEE_section_started:
            f_3.write(line)  

with open(H_in, "r") as f_2, open(INP_out, "a") as f_3:
    DistH_section_started = False
    for line in f_2:
        if Dist_start in line:
            DistH_section_started = True
        elif DistH_section_started and Dist_end in line:
            f_3.write(line)
            break
        if DistH_section_started:
            f_3.write(line)   

#Material Section
#Correction code
#define material properties for both EE and host
    #Mat1 - EE properties (always orthotropic) -> to 'correct'
with open("weave_only.inp", "r") as f_1:
    for line in f_1:
        if "type=ENGINEERING CONSTANTS" in line:
            line = f_1.readline()
            E1_e = float(line.strip().split(",")[0])
            E2_e = float(line.strip().split(",")[1])
            E3_e = float(line.strip().split(",")[2])
            v12_e = float(line.strip().split(",")[3])
            v13_e = float(line.strip().split(",")[4])
            v23_e = float(line.strip().split(",")[5])
            G12_e = float(line.strip().split(",")[6])
            G13_e = float(line.strip().split(",")[7])
            line = f_1.readline()
            G23_e = float(line.strip())
            break

#by definition the other v valus are
v32_e = (v23_e/E2_e)*E3_e; v31_e = (v13_e/E3_e)*E1_e; v21_e = (v12_e/E1_e)*E2_e;

    #Mat0 - Host properties (always isotropic) -> to use in correction, but Mat0 remains unchanged
with open("weave_only.inp", "r") as f_1:
    for line in f_1:
        if "*Material, Name=Mat0" in line:
            line = f_1.readline()
            line = f_1.readline()
            E_m = float(line.strip().split(",")[0])
            v_m = float(line.strip().split(",")[-1])
            break

#define the compliance matrices (isotropic)
C_ort_e = np.array([[1/E1_e, -v21_e/E2_e, -v31_e/E3_e, 0, 0, 0],
                 [-v12_e/E1_e, 1/E2_e, -v32_e/E3_e, 0, 0, 0],
                 [-v13_e/E1_e, -v23_e/E2_e, 1/E3_e, 0, 0, 0],
                 [0, 0, 0, 1/(2*G23_e), 0, 0],
                 [0, 0, 0, 0, 1/(2*G13_e), 0],
                 [0, 0, 0, 0, 0, 1/(2*G12_e)]]);

C_iso_m = np.array([[1/E_m, -v_m/E_m, -v_m/E_m, 0, 0, 0],
                 [-v_m/E_m, 1/E_m, -v_m/E_m, 0, 0, 0],
                 [-v_m/E_m, -v_m/E_m, 1/E_m, 0, 0, 0],
                 [0, 0, 0, (1+v_m)/E_m, 0, 0],
                 [0, 0, 0, 0, (1+v_m)/E_m, 0],
                 [0, 0, 0, 0, 0, (1+v_m)/E_m]]);

#stiffness matrix is inverse of compliance
S_ort_e = inv(C_ort_e); #original stiffness matrix for EE
S_iso_m = inv(C_iso_m); #original stiffness matrix for host matrix

#corrected stiffness matrix (s_e - s_m)
S_ort_e2 = S_ort_e - S_iso_m;

#corrected compliance matrix is inverse of corrected stiffness matrix
C_ort_e2 = inv(S_ort_e2);

#corrected youngs moduli
E1_e2 = 1/C_ort_e2[0,0];
E2_e2 = 1/C_ort_e2[1,1];
E3_e2 = 1/C_ort_e2[2,2];

#corrected poissons ratio
v12_e2 = -1*float(C_ort_e2[1,0])*E1_e2;
v13_e2 = -1*float(C_ort_e2[2,0])*E1_e2;
v23_e2 = -1*float(C_ort_e2[2,1])*E2_e2;

#corrected shear moduli
G12_e2 = 1/(C_ort_e2[3,3]);
G13_e2 = 1/(C_ort_e2[4,4]);
G23_e2 = 1/(C_ort_e2[5,5]);

#Conversion to scientfic notation
E1_e2 = "{:.3e}".format(E1_e2)
E2_e2 = "{:.3e}".format(E2_e2)
E3_e2 = "{:.3e}".format(E3_e2)

v12_e2 = "{:.3e}".format(v12_e2)
v13_e2 = "{:.3e}".format(v13_e2)
v23_e2 = "{:.3e}".format(v23_e2)

G12_e2 = "{:.3e}".format(G12_e2)
G13_e2 = "{:.3e}".format(G13_e2)
G23_e2 = "{:.3e}".format(G23_e2)

#Materials definitions (Mat0 and Mat1)
with open(INP_out,"a") as f_3:
    f_3.write("**MATERIALS\n")

    #Mat0 - always Matrix (iso)
Mat0_start = "*Material, Name=Mat0"
Mat0_end = "*Material, Name=Mat1"

with open(EE_in, "r") as f_1, open(INP_out, "a") as f_3:
    mat0_section_started = False
    for line in f_1:
        if Mat0_start in line:
            mat0_section_started = True
        elif mat0_section_started and Mat0_end in line:
            break
        if mat0_section_started:
            f_3.write(line)

    #Mat1 - always EE (ortho)
with open(INP_out, "a") as f_3:
    f_3.write("*Material, Name=Mat1\n*Elastic, type=ENGINEERING CONSTANTS\n")
    f_3.write(E1_e2 + ", " + E2_e2 + ", " + E3_e2 + ", ") #Corrected E values
    f_3.write(v12_e2 + ", " + v13_e2 + ", " + v23_e2 + ", ") #Corrected v values
    f_3.write(G12_e2 + ", " + G13_e2 + "\n" + G23_e2) #Corrected G values

#Dummy load
BC1 = "Pinned"

with open(INP_out, "a") as f_3:
    f_3.write("\n** BOUNDARY CONDITIONS\n*Boundary\nH_nset_min1, " + BC1)
    f_3.write("\n*Step, name=load\n*Static\n1., 1., 1e-05, 1.\n*Boundary\nH_nset_max1, 1, 1, 0.1\n*End Step")
