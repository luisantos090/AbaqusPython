# -*- coding: mbcs -*-
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
import os
import matplotlib.pyplot as plt
import numpy as np

# Change work directory
os.chdir(r"C:\Users\pinhosl3\Desktop\TestPython")

# VERY IMPORTANT
session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)

BaseDir = os.getcwd()
Foldername = 1

while os.path.exists(BaseDir+"/"+str(Foldername))==True:
    Foldername = Foldername + 1

os.mkdir(BaseDir+"/"+str(Foldername))
os.chdir(BaseDir+"/"+str(Foldername))



def CreateBeamModel(variables):

    top_flange_width = variables[0] #mm
    bot_flange_width = variables[1] #mm
    web_height = variables[2] #mm
    top_flange_thickness = variables[3] #mm
    bot_flange_thickness = variables[4] #mm
    web_thickness = variables[5] #mm
    Span = variables[6]

    stiffner_thickness = variables[7]
    mesh_size = variables[8]
    Load = variables[9]
    
    yield_strength = 355.0
    model_height = web_height + top_flange_thickness/2.0 + bot_flange_thickness/2.0
    
    #Create New Model
    Mdb()

    # This code creates the section through a Sketch 
    mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=200.0)
    mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0.0, 0.0), point2=(
        bot_flange_width/2.0, 0.0))
    mdb.models['Model-1'].sketches['__profile__'].HorizontalConstraint(
        addUndoState=False, entity=
        mdb.models['Model-1'].sketches['__profile__'].geometry[2])
    mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0.0, 0.0), point2=(
        -bot_flange_width/2.0, 0.0))
    mdb.models['Model-1'].sketches['__profile__'].HorizontalConstraint(
        addUndoState=False, entity=
        mdb.models['Model-1'].sketches['__profile__'].geometry[3])
    mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0.0, 0.0), point2=(
        0.0, model_height))
    mdb.models['Model-1'].sketches['__profile__'].VerticalConstraint(addUndoState=
        False, entity=mdb.models['Model-1'].sketches['__profile__'].geometry[4])
    mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0.0, model_height), point2=
        (top_flange_width/2.0, model_height))
    mdb.models['Model-1'].sketches['__profile__'].HorizontalConstraint(
        addUndoState=False, entity=
        mdb.models['Model-1'].sketches['__profile__'].geometry[5])
    mdb.models['Model-1'].sketches['__profile__'].PerpendicularConstraint(
        addUndoState=False, entity1=
        mdb.models['Model-1'].sketches['__profile__'].geometry[4], entity2=
        mdb.models['Model-1'].sketches['__profile__'].geometry[5])
    mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0.0, model_height), point2=
        (-top_flange_width/2.0, model_height))
    mdb.models['Model-1'].sketches['__profile__'].HorizontalConstraint(
        addUndoState=False, entity=
        mdb.models['Model-1'].sketches['__profile__'].geometry[6])
    mdb.models['Model-1'].Part(dimensionality=THREE_D, name='Part-1', type=
        DEFORMABLE_BODY)
    mdb.models['Model-1'].parts['Part-1'].BaseShellExtrude(depth=Span, sketch=
        mdb.models['Model-1'].sketches['__profile__'])
    del mdb.models['Model-1'].sketches['__profile__']


    #This creates the stiffner geometry
    mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=200.0)
    mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0.0, 0.0), point2=(
        bot_flange_width/2.0, 0.0))
    mdb.models['Model-1'].sketches['__profile__'].HorizontalConstraint(
        addUndoState=False, entity=
        mdb.models['Model-1'].sketches['__profile__'].geometry[2])
    mdb.models['Model-1'].sketches['__profile__'].Line(point1=(bot_flange_width/2.0, 0.0), point2=(
        top_flange_width/2.0, model_height))
    mdb.models['Model-1'].sketches['__profile__'].Line(point1=(top_flange_width/2.0, model_height), point2=
        (-top_flange_width/2.0, model_height))
    mdb.models['Model-1'].sketches['__profile__'].HorizontalConstraint(
        addUndoState=False, entity=
        mdb.models['Model-1'].sketches['__profile__'].geometry[4])
    mdb.models['Model-1'].sketches['__profile__'].Line(point1=(-top_flange_width/2.0, model_height), 
        point2=(-bot_flange_width/2.0, 0.0))
    mdb.models['Model-1'].sketches['__profile__'].Line(point1=(-bot_flange_width/2.0, 0.0), point2=
        (0.0, 0.0))
    mdb.models['Model-1'].sketches['__profile__'].HorizontalConstraint(
        addUndoState=False, entity=
        mdb.models['Model-1'].sketches['__profile__'].geometry[6])
    mdb.models['Model-1'].Part(dimensionality=THREE_D, name='Part-2', type=
        DEFORMABLE_BODY)
    mdb.models['Model-1'].parts['Part-2'].BaseShell(sketch=
        mdb.models['Model-1'].sketches['__profile__'])
    del mdb.models['Model-1'].sketches['__profile__']


    # This creates the material
    mdb.models['Model-1'].Material(name='Material-1')
    mdb.models['Model-1'].materials['Material-1'].Elastic(table=((210000.0, 0.3), 
        ))
    mdb.models['Model-1'].materials['Material-1'].Plastic(table=((yield_strength, 0.0), ))


    # This creates the 3 sections with different thicknesses 
    mdb.models['Model-1'].HomogeneousShellSection(idealization=NO_IDEALIZATION, 
        integrationRule=SIMPSON, material='Material-1', name='Section-1', 
        nodalThicknessField='', numIntPts=5, poissonDefinition=DEFAULT, 
        preIntegrate=OFF, temperature=GRADIENT, thickness=top_flange_thickness, thicknessField='', 
        thicknessModulus=None, thicknessType=UNIFORM, useDensity=OFF)
    mdb.models['Model-1'].HomogeneousShellSection(idealization=NO_IDEALIZATION, 
        integrationRule=SIMPSON, material='Material-1', name='Section-2', 
        nodalThicknessField='', numIntPts=5, poissonDefinition=DEFAULT, 
        preIntegrate=OFF, temperature=GRADIENT, thickness=web_thickness, thicknessField='', 
        thicknessModulus=None, thicknessType=UNIFORM, useDensity=OFF)
    mdb.models['Model-1'].HomogeneousShellSection(idealization=NO_IDEALIZATION, 
        integrationRule=SIMPSON, material='Material-1', name='Section-3', 
        nodalThicknessField='', numIntPts=5, poissonDefinition=DEFAULT, 
        preIntegrate=OFF, temperature=GRADIENT, thickness=bot_flange_thickness, thicknessField='', 
        thicknessModulus=None, thicknessType=UNIFORM, useDensity=OFF)
    mdb.models['Model-1'].HomogeneousShellSection(idealization=NO_IDEALIZATION, 
        integrationRule=SIMPSON, material='Material-1', name='Section-4', 
        nodalThicknessField='', numIntPts=5, poissonDefinition=DEFAULT, 
        preIntegrate=OFF, temperature=GRADIENT, thickness=stiffner_thickness, thicknessField='', 
        thicknessModulus=None, thicknessType=UNIFORM, useDensity=OFF)
        
    # This assigns the sections to each plates
    # This was done in the following order: Top flange, web, bottom flange and stiffner
    mdb.models['Model-1'].parts['Part-1'].SectionAssignment(offset=0.0, 
        offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
        faces=mdb.models['Model-1'].parts['Part-1'].faces.getSequenceFromMask(
        mask=('[#11 ]', ), )), sectionName='Section-1', thicknessAssignment=
        FROM_SECTION)
    mdb.models['Model-1'].parts['Part-1'].SectionAssignment(offset=0.0, 
        offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
        faces=mdb.models['Model-1'].parts['Part-1'].faces.getSequenceFromMask(
        mask=('[#2 ]', ), )), sectionName='Section-2', thicknessAssignment=
        FROM_SECTION)
    mdb.models['Model-1'].parts['Part-1'].SectionAssignment(offset=0.0, 
        offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
        faces=mdb.models['Model-1'].parts['Part-1'].faces.getSequenceFromMask(
        mask=('[#c ]', ), )), sectionName='Section-3', thicknessAssignment=
        FROM_SECTION)
    mdb.models['Model-1'].parts['Part-2'].SectionAssignment(offset=0.0, 
        offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
        faces=mdb.models['Model-1'].parts['Part-2'].faces.getSequenceFromMask(
        mask=('[#1 ]', ), )), sectionName='Section-4', thicknessAssignment=
        FROM_SECTION)
        
    # This starts the assemble
    mdb.models['Model-1'].rootAssembly.DatumCsysByDefault(CARTESIAN)
    mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='Part-1-1', 
        part=mdb.models['Model-1'].parts['Part-1'])
    mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='Part-2-1', 
        part=mdb.models['Model-1'].parts['Part-2'])
        
    # This translates the stiffner to mid-span    
    mdb.models['Model-1'].rootAssembly.translate(instanceList=('Part-2-1', ), 
        vector=(0.0, 0.0, Span/2.0))

    # This merges the two parts into one
    mdb.models['Model-1'].rootAssembly.InstanceFromBooleanMerge(domain=GEOMETRY, 
        instances=(mdb.models['Model-1'].rootAssembly.instances['Part-1-1'], 
        mdb.models['Model-1'].rootAssembly.instances['Part-2-1']), name='Part-3', 
        originalInstances=SUPPRESS)
        
    # This creates the step
    mdb.models['Model-1'].StaticStep(initialInc=0.01, maxInc=0.01, name='Step-1', 
        nlgeom=ON, previous='Initial')
        
        
    # This applies the load
    mdb.models['Model-1'].ConcentratedForce(cf2=-Load, createStepName='Step-1', 
        distributionType=UNIFORM, field='', localCsys=None, name='Load-1', region=
        Region(
        vertices=mdb.models['Model-1'].rootAssembly.instances['Part-3-1'].vertices.findAt(
        ((0.0, model_height, Span/2.0), ), )))
        
    # This is the boundary conditions

    mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
        distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
        'BC-1', region=Region(
        edges=mdb.models['Model-1'].rootAssembly.instances['Part-3-1'].edges.findAt(
        ((-bot_flange_width/4.0, 0.0, Span), ), ((bot_flange_width/4.0, 0.0, Span), ), )), u1=0.0, u2=0.0, u3=
        0.0, ur1=UNSET, ur2=UNSET, ur3=UNSET)
    mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
        distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
        'BC-2', region=Region(
        edges=mdb.models['Model-1'].rootAssembly.instances['Part-3-1'].edges.findAt(
        ((bot_flange_width/4.0, 0.0, 0.0), ), ((-bot_flange_width/4.0, 0.0, 0.0), ), )), u1=0.0, u2=0.0, u3=UNSET, 
        ur1=UNSET, ur2=UNSET, ur3=UNSET)

        
    # This meshes the assembly
    mdb.models['Model-1'].parts['Part-3'].seedPart(deviationFactor=0.1, 
        minSizeFactor=0.1, size=mesh_size)
    mdb.models['Model-1'].parts['Part-3'].generateMesh()

        
        
    # This creates the job!
    mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
        explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
        memory=90, memoryUnits=PERCENTAGE, model='Model-1', modelPrint=OFF, 
        multiprocessingMode=DEFAULT, name='Job-1', nodalOutputPrecision=SINGLE, 
        numCpus=1, numGPUs=0, queue=None, resultsFormat=ODB, scratch='', type=
        ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)

    # This submits the job
    mdb.jobs['Job-1'].submit(consistencyChecking=OFF)


    # TO WAIT FOR JOB COMPLETION
    mdb.jobs['Job-1'].waitForCompletion()
    print("SS I-beam Model finished running")
    

def PostProcessing():
    CurrentDir = os.getcwd()
    odb = session.openOdb(CurrentDir + '/Job-1.odb')
    NrOfSteps = len(odb.steps['Step-1'].frames)
    
    displacements = []
    
    for i in range(NrOfSteps):
        central_disp = odb.steps['Step-1'].frames[i].fieldOutputs['U'].values[0].data[1]*-1
        displacements.append(central_disp)
        
    Forces = []
    
    for i in range(NrOfSteps):
        applied_force = odb.steps['Step-1'].frames[i].fieldOutputs['CF'].values[1].data[1]*-1
        Forces.append(applied_force)
 
    fig, ax = plt.subplots()
    ax.plot(displacements, Forces, color='r', label='U2')
    plt.legend()
    ax.set(xlabel='Displacements [mm]', ylabel='Applied Load [kN]',
           title='Force Displacement Curve')
    ax.grid()

    fig.savefig("MAX_DISPLACMENT.png")
    plt.close(fig)

"""
top_flange_width = variables[0] #mm
bot_flange_width = variables[1] #mm
web_height = variables[2] #mm
top_flange_thickness = variables[3] #mm
bot_flange_thickness = variables[4] #mm
web_thickness = variables[5] #mm
Span = variables[6]

stiffner_thickness = variables[7]
mesh_size = variables[8]
Load = variables[9]
"""
   
models = []

models.append([180.0, 200.0, 190.0, 8.0, 9.0, 10.0, 9000.0, 11.0, 100.0, 5000.0])
models.append([180.0, 200.0, 190.0, 8.0, 9.0, 10.0, 8000.0, 11.0, 100.0, 5000.0])
models.append([180.0, 200.0, 190.0, 8.0, 9.0, 10.0, 7000.0, 11.0, 100.0, 5000.0])
models.append([180.0, 200.0, 190.0, 8.0, 9.0, 10.0, 6000.0, 11.0, 100.0, 5000.0])
models.append([180.0, 200.0, 190.0, 8.0, 9.0, 10.0, 5000.0, 11.0, 100.0, 5000.0])
models.append([180.0, 200.0, 190.0, 8.0, 9.0, 10.0, 4000.0, 11.0, 100.0, 5000.0])

for i in models:
    ModelName = str(i[0])
    for j in i[1:]:
        ModelName = ModelName + "_" + str(j)
    NameFolder = BaseDir+ "/" +str(Foldername) + "/" + ModelName 
    print(NameFolder)
    os.mkdir(NameFolder)
    os.chdir(NameFolder)
    CreateBeamModel(i)
    PostProcessing()