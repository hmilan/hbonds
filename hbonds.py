#! /usr/bin/python
# probably works for 3 only ??

import sys,time

if sys.argv[1:]: inpf = open(sys.argv[1],mode='rt',encoding='utf-8')
else: sys.stderr.write('needs input file name?\n') ; exit()

# OLD:
#if sys.argv[2:]: outf = open(sys.argv[2],mode='wt',encoding='utf-8')
#else:            outf = sys.stdout

# these 2 lines maybe *CHANGE* -d
# we chose the atoms that are selected on the "coor hbond ..." command
# in our case i and j are the same!
# additionally if i==j then this column must be zero!
# generally these lists can be different, both names and lengths
atom_type_i=['HNF','OH2','OF','O3','O11','O12','O13','O14']
atom_type_j=[x for x in atom_type_i] 
# It works also without the following line,
# but don't clutter the folder with too many files:
# water is not in the first selection
atom_type_i.remove('OH2')

#OLD:
#atom_type_i=['HNF']
#atom_type_j=['OH2','OF','O3','O11','O12','O13','O14']

i_types=len(atom_type_i) #; print('n_types=',n_types)
j_types=len(atom_type_j) #; print('n_types=',n_types)
# results=(n_types+1)*[[]] #- WEIRD this definition is very bad ???
# WEIRD: next 2 lines work OK
#results=[]
#for i in range(n_types+1): results.append([])
# OR better:
results=[[[] for j in range(j_types+1)] for i in range(i_types)]
# +1 is for total over j atoms!
# results[][][] - first list is for atom_i second is for atom_j and last goes through the timeline

# *CHANGE* -able: t_step (usually related to frame saving time)
# this could be read from output!!
t_step=5.0  # unit=ps IMPORTANT!!! Must be correct

count=i_types*[0] ; life_time=i_types*[0.0]
curr_time=i_types*[0.0] ; prev_time=i_types*[0.0]
curr_index=i_types*[0]

time_series_flag=False
occupancy_flag=False

for line in inpf:
    if line.find('(Bridge)') > 1:
        time_series_flag=True
        #outf.write(line)
    if line.find('<occupancy>') > 1:
        occupancy_flag=True
        time_series_flag=False
        #outf.write(line)
    if time_series_flag:
        # now we are collecting data for trajectory info
        if len(line) > 1:
            # *CHANGE* - able - just to ignore non-relevant lines
            if line.find('MEMB'):
#                if line.split()[3].find(atom_type_i[0]) == 0 :
                if line.split()[3].strip() in atom_type_i :
                    iatom=atom_type_i.index(line.split()[3].strip())
                    # some of these [iatom] lists are maybe not needed ??
                    # but they don't hurt either
                    count[iatom] += 1
                    curr_time[iatom]=float(line.split()[10])
                    life_time[iatom]=float(line.split()[9])
                    # make space for next time step
                    if curr_time[iatom] > prev_time[iatom]:
                        prev_time[iatom]=curr_time[iatom]
                        for j_type in range(j_types+1):
                            results[iatom][j_type].append(float(0.0))

                    if line.split()[8] in atom_type_j:
                        atom=atom_type_j.index(line.split()[8])
                    else: continue  # we are not interested in this atom...
                    curr_index[iatom]=len(results[iatom][atom])
                    for tt in range(int(life_time[iatom]/t_step),0,-1):
                        results[iatom][atom][curr_index[iatom]-tt] += 1

    if occupancy_flag : break

                #if count > 300 : break

#print(results)

# last column is the sum of all previous columns
#OLD:results[j_types]=[sum([a[i] for a in results[:-1]]) for i in range(len(results[0]))]
for ix in range(i_types):
#    print('times: ',len(results[ix][0]))
    for i in range(len(results[ix][0])):
        s=0.0
        for k in range(len(results[ix])): s += results[ix][k][i]
        results[ix][j_types][i]=s  #sum([a for a in results[ix][:-1][i]])
#        print(results[ix][j_types][i])

# write all files for atom_type_i[] atoms
for atname in atom_type_i:
    outf = open(sys.argv[1]+'-'+atname+'.dat',mode='wt',encoding='utf-8')
    iatom=atom_type_i.index(atname)
    # *CHANGE* Header for plot file
    outf.write('# count = {0:d}'.format(count[iatom])+" for atom: "+atname+"\n")
    outf.write("#   time       ");
    for i in range(len(atom_type_j)): outf.write("{0:6s} ".format(atom_type_j[i]))
    outf.write("total\n");

    # time series for the counts...
    # *CHANGE* start time, units, etc...
    # normalize occupancy (eg 288 for 12x12 lipids)
    t_start=2*t_step
    t_units=1.0 # ps, for ns put 0.001
    t_units=0.001 # use ns here!
    t_occ=1.0 # no normalization
    for i in range(len(results[iatom][0])):
        outf.write('{0:12.3f} '.format((i*t_step+t_start)*t_units))
        for j in range(len(results[iatom])):
            outf.write("{0:7.0f}".format(results[iatom][j][i]*t_occ))
        outf.write('\n')

    outf.close()
