from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
import math
import numpy as np
import matplotlib.pyplot as plt
import os
import random
import shutil
import pickle


# Change work directory
path = r"C:\Users\pinhosl3\Documents\GA_video"
os.chdir(path)


# VERY IMPORTANT
session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)


def PrintToScreen(Message):
    print >> sys.__stdout__, '%s' % Message
    


length = 6000.0
width = 5000.0


No_columns = 4

x_range = [200.0, length - 200.0, 8]
y_range = [200.0, width - 200.0, 7]

for i in [x_range, y_range]:
    step = (i[1]-i[0])/(2**i[2]-1)
    i.append(step)

Population_size = 100
Generations = 25
Elitism_size = 2
Tournament_size = 2 

variables = []
for i in range(No_columns):
    variables.append(x_range)
    variables.append(y_range)
    
PrintToScreen('New Run[write "N"]? Or Continue Previous Run[write "C"]?')
GA_status = input('New Run[write "N"]? Or Continue Previous Run[write "C"]?')
if str(GA_status) == "N":
    GA_status = "New"
elif str(GA_status) == "C":
    GA_status = "Continue"


if GA_status == "Continue":
    Round = 1
    while os.path.isdir("%s/Round%s" % (path, Round + 1)):
        Round = Round + 1
        PrintToScreen(Round)
    os.chdir("%s/Round%s" % (path, Round))
    for dirs in os.listdir("%s/Round%s" % (path, Round)):
        if dirs.startswith("[") and os.path.isfile("%s/Round%s/%s/FitnessData.txt" % (path, Round, dirs)) == False:
            shutil.rmtree("%s/Round%s/%s" % (path, Round, dirs))
            break
    

def CreateInitialPop(variables, Population_size):
    count = 0 
    InitPop = []
    Chromossome_Lenght = 0
    for i in variables:
        Chromossome_Lenght = Chromossome_Lenght + i[2]
    while count < Population_size:
        Chromossome = ''
        for i in range(Chromossome_Lenght):
            Chromossome = Chromossome + str(random.choice([0,1]))        
        if Chromossome not in InitPop:
            InitPop.append([Chromossome])
            count = count + 1
        else:
            print("------REPEATED-------")
    
    return InitPop
    
    
def ChromossomeToData(Chromossome, variables):
    var_data = []
    start = 0
    col_data = []
    for i in variables:
        var_binary = Chromossome[0][start:start+i[2]]
        var_decimal = i[0] + int(str(var_binary), 2)*i[3]
        col_data.append(round(var_decimal,1))
        start = start + i[2]
        if len(col_data) == 2:
            var_data.append(col_data)
            col_data = []
    return var_data


def CalculateDisplacement(Columns_coordinates):
    # Create a new model
    Mdb()

    # This creates a rectangular slab
    mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=200.0)
    mdb.models['Model-1'].sketches['__profile__'].rectangle(point1=(0.0, 0.0), 
        point2=(length, width))
    mdb.models['Model-1'].Part(dimensionality=THREE_D, name='Part-1', type=
        DEFORMABLE_BODY)
    mdb.models['Model-1'].parts['Part-1'].BaseShell(sketch=
        mdb.models['Model-1'].sketches['__profile__'])
    del mdb.models['Model-1'].sketches['__profile__']

    # This creates the material - Assumes linear elastic behaviour 
    mdb.models['Model-1'].Material(name='Material-1')
    mdb.models['Model-1'].materials['Material-1'].Elastic(table=((30000.0, 0.2), ))

    # Creates section
    mdb.models['Model-1'].HomogeneousShellSection(idealization=NO_IDEALIZATION, 
        integrationRule=SIMPSON, material='Material-1', name='Section-1', 
        nodalThicknessField='', numIntPts=5, poissonDefinition=DEFAULT, 
        preIntegrate=OFF, temperature=GRADIENT, thickness=300.0, thicknessField='', 
        thicknessModulus=None, thicknessType=UNIFORM, useDensity=OFF)
    mdb.models['Model-1'].parts['Part-1'].SectionAssignment(offset=0.0, 
        offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
        faces=mdb.models['Model-1'].parts['Part-1'].faces.findAt(((1.0, 
        1.0, 0.0), (0.0, 0.0, 1.0)), )), sectionName='Section-1', 
        thicknessAssignment=FROM_SECTION)
        
    # This creates the assembly
    mdb.models['Model-1'].rootAssembly.DatumCsysByDefault(CARTESIAN)
    mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='Part-1-1', 
        part=mdb.models['Model-1'].parts['Part-1'])
    mdb.models['Model-1'].rootAssembly.makeIndependent(instances=(
        mdb.models['Model-1'].rootAssembly.instances['Part-1-1'], ))
        
    # This sketches the supports
    mdb.models['Model-1'].ConstrainedSketch(gridSpacing=390.51, name='__profile__', 
        sheetSize=15620.49, transform=
        mdb.models['Model-1'].rootAssembly.MakeSketchTransform(
        sketchPlane=mdb.models['Model-1'].rootAssembly.instances['Part-1-1'].faces.findAt(
        (1.0, 1.0, 0.0), (0.0, 0.0, 1.0)), sketchPlaneSide=SIDE1, 
        sketchUpEdge=mdb.models['Model-1'].rootAssembly.instances['Part-1-1'].edges.findAt(
        (length, 1.0, 0.0), ), sketchOrientation=RIGHT, origin=(0.0, 0.0, 0.0)))
    mdb.models['Model-1'].rootAssembly.projectReferencesOntoSketch(filter=
        COPLANAR_EDGES, sketch=mdb.models['Model-1'].sketches['__profile__'])


    for col_coord in Columns_coordinates:
        mdb.models['Model-1'].sketches['__profile__'].rectangle(point1=(col_coord[0]-150.0, col_coord[1]-150.0), point2=(col_coord[0]+150.0, col_coord[1]+150.0))


    mdb.models['Model-1'].rootAssembly.PartitionFaceBySketch(faces=
        mdb.models['Model-1'].rootAssembly.instances['Part-1-1'].faces.findAt(((
        1.0, 1.0, 0.0), )), sketch=
        mdb.models['Model-1'].sketches['__profile__'], sketchUpEdge=
        mdb.models['Model-1'].rootAssembly.instances['Part-1-1'].edges.findAt((
        length, 1.0, 0.0), ))
    del mdb.models['Model-1'].sketches['__profile__']


    # This creates the step - 1 iteration
    mdb.models['Model-1'].StaticStep(initialInc=1.0, maxInc=1.0, name='Step-1', 
        previous='Initial')


    Columns = mdb.models['Model-1'].rootAssembly.instances['Part-1-1'].faces.findAt(((Columns_coordinates[0][0], Columns_coordinates[0][1], 0.0), (0.0, 0.0, 1.0)))
    for col_coord in Columns_coordinates[1:]:
        Columns = Columns + mdb.models['Model-1'].rootAssembly.instances['Part-1-1'].faces.findAt(((col_coord[0], col_coord[1], 0.0), (0.0, 0.0, 1.0)))


    Slab = mdb.models['Model-1'].rootAssembly.instances['Part-1-1'].faces.findAt(((1.0, 1.0, 0.0), (0.0, 0.0, 1.0)))

    # This creates the Pressure
    mdb.models['Model-1'].Pressure(amplitude=UNSET, createStepName='Step-1', 
        distributionType=UNIFORM, field='', magnitude=1.0, name='Load-1', region=
        Region(side1Faces=Columns+Slab))
        
    mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
        distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
        'BC-1', region=Region(faces=Columns), u1=
        0.0, u2=0.0, u3=0.0, ur1=UNSET, ur2=UNSET, ur3=UNSET)
        
    # Mesh creation
    mdb.models['Model-1'].rootAssembly.seedPartInstance(deviationFactor=0.1, 
        minSizeFactor=0.1, regions=(
        mdb.models['Model-1'].rootAssembly.instances['Part-1-1'], ), size=100.0)
    mdb.models['Model-1'].rootAssembly.generateMesh(regions=(
        mdb.models['Model-1'].rootAssembly.instances['Part-1-1'], ))
        
    # This creates the Job
    mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
        explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
        memory=90, memoryUnits=PERCENTAGE, model='Model-1', modelPrint=OFF, 
        multiprocessingMode=DEFAULT, name='Job-1', nodalOutputPrecision=SINGLE, 
        numCpus=1, numGPUs=0, queue=None, resultsFormat=ODB, scratch='', type=
        ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)
    
    # SUBMIT THE JOB
    mdb.jobs['Job-1'].submit(consistencyChecking=OFF)
    # TO WAIT FOR JOB COMPLETION
    mdb.jobs['Job-1'].waitForCompletion()
    print("Job-1 finished running")
    
    CurrentFolder = os.getcwd()
    
    # This opens the odb - Output Database
    odb = session.openOdb(str(CurrentFolder)+'/'+'Job-1.odb')
    
    # This get the displacements
    # frame = 0 is the initial step - frame = 1 is the first fully calculated step
    Disps = odb.steps['Step-1'].frames[1].fieldOutputs['U'].getSubset(position = NODAL).bulkDataBlocks[0].data
    MaxDisps = np.max(np.abs(Disps),axis=0)
    MaxU3 = MaxDisps[2]
    
    odb.close()
    print(MaxU3)
    
    if os.path.exists("%s/%s" % (str(CurrentFolder), "Job-1.simdir")):
        shutil.rmtree("%s/%s" % (str(CurrentFolder), "Job-1.simdir"))
    
    for fname in os.listdir(CurrentFolder):
        if fname.startswith("Job-1") and not fname.endswith(".inp"):
            os.remove(os.path.join(CurrentFolder, fname))
    
    return MaxU3


def Tournament(OldPop, PlayersNumber):
    ww = np.random.randint(0,len(OldPop))
    winner = OldPop[ww][0]
    counter = 0
    while counter < PlayersNumber - 1:
        cc = np.random.randint(0,len(OldPop))
        challenger = OldPop[cc][0]
        #Remember to verify that the Fitness is the last value
        # Check later if we want to minimise or maximise the fitness - I think this is correct already
        if OldPop[ww][-1] < OldPop[cc][-1]:
            continue
        else:
            ww = cc
            winner = challenger
        counter = counter + 1
    return winner


def NewPopulation(OldPop, select, Tournament_size, Current_Round, Generations):
    Newpop = []
    counter = select
    OldPop_sorted = sorted(OldPop,key=lambda x: x[-1])
    # This does elitism selection
    for i in range(select):
        Newpop.append([OldPop_sorted[i][0]])
    # This does tournament selection before crossover
    while counter < len(OldPop):
        father = Tournament(OldPop, Tournament_size)
        mother = Tournament(OldPop, Tournament_size)
        # while mother == father:
            # mother = Tournament(OldPop, fitness, 3)
        child1 = ''
        child2 = ''
        # This does crossover of genes: 70% chance of occuring
        if np.random.uniform(0, 1) <= 0.7:
            for i in range(len(OldPop[0][0])):    
                child1 = child1 + str(np.random.choice([father[i],mother[i]]))
                child2 = child2 + str(np.random.choice([father[i],mother[i]]))
        else:
            child1 = father
            child2 = mother
        counter = counter + 2
        # This does mutation of genes
        # child1 = MutateChromossome(child1, 0.1)
        # child2 = MutateChromossome(child2, 0.1)
        child1 = MutateChromossome(child1, 0.1+0.8*(Current_Round/float(Generations)))
        child2 = MutateChromossome(child2, 0.1+0.8*(Current_Round/float(Generations)))
        Newpop.append([child1])
        Newpop.append([child2])
    return Newpop


def MutateChromossome(Chromossome, Probability):
    if np.random.uniform(0, 1) <= Probability:
        Random_Gene = np.random.randint(0,len(Chromossome))
        if Chromossome[Random_Gene] == '0':
            Chromossome = Chromossome[0:Random_Gene] + '1' + Chromossome[Random_Gene+1:]
        else:
            Chromossome = Chromossome[0:Random_Gene] + '0' + Chromossome[Random_Gene+1:]
    return Chromossome


def CreateFolder(name):
    if not os.path.exists("%s" % name):
        os.makedirs("%s" % name)
    os.chdir("%s" % name)


def GeneticAlgorithm(variables, Generations, Population_size, Elitism_size, Tournament_size, Round):
    if Round == 1:
        Process = []
        os.chdir("%s/Round%s" % (path,Round))
        Population_binary = pickle.load(open("Population_binary_list", 'rb'))
        SaveAnalysisbinary = pickle.load(open("SaveAnalysisbinary", 'rb'))
    elif Round == 0:
        Process = []
        SaveAnalysisbinary = dict()
        Population_binary = CreateInitialPop(variables, Population_size)
        Round = 1
    else:
        os.chdir(path)
        Process = pickle.load(open("Process_list", 'rb'))
        os.chdir("%s/Round%s" % (path,Round))
        Population_binary = pickle.load(open("Population_binary_list", 'rb'))
        SaveAnalysisbinary = pickle.load(open("SaveAnalysisbinary", 'rb'))
        
    for i in range(Generations-Round+1):
        os.chdir(path)
        CreateFolder("Round%s" % (i+Round))
        PrintToScreen('-------ROUND %s-------' % (i+Round))
        os.chdir("%s/Round%s" % (path,i+Round))
        pickle.dump(Population_binary, open("Population_binary_list", 'wb'))
        if len(Population_binary[0])==1:
            for j in range(Population_size):
                Population_binary[j].append(ChromossomeToData(Population_binary[j], variables))
        for Chromossome in Population_binary:
            if str(Chromossome[0]) in SaveAnalysisbinary:
                PrintToScreen("Saving analysis")
                #print("Retrieving results from previous analysis")
                Fitness = SaveAnalysisbinary[str(Chromossome[0])][-1]
                Chromossome.append(Fitness)
                os.chdir("%s/Round%s" % (path,i+Round))
            else:
                #print("New binary")
                PrintToScreen("New binary")
                CreateFolder("%s" % (Chromossome[1]))
                os.chdir("%s/Round%s/%s" % (path,i+Round,Chromossome[1]))
                Fitness = CalculateDisplacement(Chromossome[1])
                SaveAnalysisbinary[str(Chromossome[0])] = [Chromossome[1], (i+Round), Fitness]
                Chromossome.append(Fitness)
                os.chdir("%s/Round%s" % (path,i+Round))
                pickle.dump(SaveAnalysisbinary, open("SaveAnalysisbinary", 'wb'))
        os.chdir(path)
        MinFitness = 1.0e8 
        SumFitness = 0
        for Chromossome in Population_binary:
            SumFitness = SumFitness + Chromossome[-1]
            MinFitness = min(Chromossome[-1], MinFitness)
        Process.append([i+Round, MinFitness, SumFitness/Population_size])
        PrintData("%s/Round%s" % (path,i+Round), 'Process', "Round\tMinFitness\tAverageFitness", Process)
        PrintData("%s/Round%s" % (path,i+Round), 'Population_binary', "GeneticCode\tVariables\tFitness", Population_binary)
        PrintFinalData("%s/Round%s" % (path,i+Round), 'AllIndividuals', "GeneticCode\tVariables\tRound\tFitness", SaveAnalysisbinary)
        pickle.dump(Process, open("Process_list", 'wb'))
        if i+Round == Generations-1:
            continue
        else:
            Population_binary = NewPopulation(Population_binary, Elitism_size, Tournament_size, i+1, Generations)
    
    PrintData("%s/Round%s" % (path,i+Round), 'Process', "Round\tMinFitness\tAverageFitness", Process)
    PrintFinalData("%s/Round%s" % (path,i+Round), 'AllIndividuals', "GeneticCode\tVariables\tRound\tFitness", SaveAnalysisbinary)
    #print('GA is completed')
    PrintToScreen('GA is completed')


def DeleteRounds(path):
    os.chdir(path)
    for dirs in os.listdir(path):
        if dirs.startswith("Round"):
            shutil.rmtree("%s/%s" % (path,dirs))


def PrintData(Folder, Filename, Title, Data):
    #This will depend on the application in mind
    opFile = Folder+"/"+Filename+'.txt'
    
    try:
        opFileU = open(opFile,'w')
        opFileU.write("%10s\n" % Title)
    except IOError:
        print('cannot open', opFile)
        exit(0)

    for line in Data:
        for item in line:
            opFileU.write(str(item))
            opFileU.write("\t")
        opFileU.write("\n")
    opFileU.close()


def PrintFinalData(Folder, Filename, Title, Data):
    #This will depend on the application in mind
    opFile = Folder+"/"+Filename+'.txt'
    
    try:
        opFileU = open(opFile,'w')
        opFileU.write("%10s\n" % Title)
    except IOError:
        print('cannot open', opFile)
        exit(0)

    GenCodes = list(Data.keys())
    Fitnesses = list(Data.values())
    for i in range(len(Data)):
        opFileU.write(GenCodes[i])
        opFileU.write("\t")
        opFileU.write(str(Fitnesses[i]))
        # opFileU.write(Chromossome[str(Chromossome[0])])
        opFileU.write("\n")
    opFileU.close()

    
if GA_status == "New":
    DeleteRounds("%s" % path)
    Round = 0
    GeneticAlgorithm(variables, Generations, Population_size, Elitism_size, Tournament_size, Round)
    
elif GA_status == "Continue":
    GeneticAlgorithm(variables, Generations, Population_size, Elitism_size, Tournament_size, Round)



     
# Columns_coordinates = []
# Columns_coordinates.append([500.0,600.0])
# Columns_coordinates.append([700.0,4200.0])
# Columns_coordinates.append([4500.0,900.0])
# Columns_coordinates.append([5200.0,3900.0])

# MaxU3 = CalculateDisplacement(Columns_coordinates)
# print(MaxU3)


