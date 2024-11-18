filedisp = 'D:\iFEM\\temp\\ODB2VTK-main\\temp\Job-2\\'
S_iFEMfile = filedisp + 'file2_Ndsp\\SBT_1.txt'
with open(S_iFEMfile, 'r') as stressfile:
    U_ERR = stressfile.readlines()
    print(U_ERR)
    for line in U_ERR:
        print(line)
        result = line.split("\t")
        print(result)

a = [1,2,3]
a[1] = 5
print(a)