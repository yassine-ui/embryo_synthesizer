import sys
import vtk
import numpy as np
import math
import os
import pathlib
import subprocess
from vtk.util.misc import vtkGetDataRoot
from vtk.util.colors import peacock, tomato
def main():
    PACKINGS_PATH = "./configs/packings/"
    COUNTER = 0
    for f in os.listdir(PACKINGS_PATH):
        if "opt" in f:
            try:
                # Détermination du packing d'ellipsides spécifiés dans un fichier externe (.opt)
                fileIn = os.path.join(PACKINGS_PATH, f)
                args=sys.argv
                if len(args)>1:
                    fileIn=args[1]
                if len(args)>2:
                    pz=np.float(args[2])
                else:
                    pz=0.5
                subprocess.call(["./packEllipsoids",fileIn])
                COUNTER += 1
                # Lecture de la solution de packing (algo packingellipsoids [Lobato])
                (RSphere,NbEllipsoides,Ellipsoides)=ReadResPacking('solution.txt')

                # Mise en équation des ellipsoides et de la sphère englobante

                #(A,S)=FromGeomToCartForEllipsoid(NbEllipsoides,Ellipsoides,RSphere)
                #print(A)
                #lnume=np.array(range(NbEllipsoides))
                volumes=GetVolumes(NbEllipsoides,Ellipsoides)
                #print(volumes)

                # Affichage es ellipsoides et de la sphère englobante
                renWin=PlotEllipsoideSimple(RSphere, NbEllipsoides, Ellipsoides, pz)

                #Création d'une image présentant les contours de l'intersection des ellipsoides avec un plan
                CreateImageFromScene(renWin,"image_"+str(COUNTER)+".png")
            except Exception as e:
                print(e)

def ReadResPacking(filename):
# lecture de la sortie du programme de packing dans une sphère

    with open(filename, 'r') as infile:
        data = infile.read()
    my_list = data.splitlines()
    l = 0
    print(my_list[l])
    RSphere=np.float(my_list[0])
    NbEllipsoides=np.int(my_list[1])
    Ellipsoides= [float(x) for x in ' '.join(my_list[2:]).split()]
    Ellipsoides= np.reshape(Ellipsoides,(NbEllipsoides,9))
    infile.close()
    return(RSphere,NbEllipsoides,Ellipsoides)

def FromGeomToCartForEllipsoid(NbEllipsoides,Ellipsoides,RSphere):
    #LDemiAxes,Centre,Angles
    # Passage des demi longueurs d'axes, du centre et de l'angle de rotation, phi (resp. theta, psi)  autor de x (resp., y, Z)
    # à une équation en cartésiennes : a0x2+a1y2+a2z2+a3xy+a4yz+a5xz+a6x+a7y+a8z+a9=0
    import numpy as np
    a=Ellipsoides[:,0]
    b=Ellipsoides[:,1]
    c=Ellipsoides[:,2]
    centres=Ellipsoides[:,3:6]
    phi=Ellipsoides[:,6]
    theta=Ellipsoides[:,7]
    psi=Ellipsoides[:,8]
    cphi=np.cos(phi)
    sphi=np.sin(phi)
    ctheta=np.cos(theta)
    stheta=np.sin(theta)
    cpsi=np.cos(psi)
    spsi=np.sin(psi)
    Q1=[ctheta*cpsi, sphi*stheta*cpsi-cphi*spsi, sphi*spsi+cphi*stheta*cpsi]
    Q2=[ctheta*spsi, cphi*cpsi+sphi*stheta*spsi, cphi*stheta*spsi- sphi*cpsi]
    Q3=[-stheta, sphi*ctheta, cphi*ctheta]
    Q=[Q1,Q2,Q3]
    Z=np.zeros(NbEllipsoides)
    A=np.zeros([NbEllipsoides,10])
    PI=[[1/a,Z,Z],[Z,1/b,Z],[Z,Z,1/c]]
    for k in range(NbEllipsoides):
        pe=[[PI[0][0][k],PI[0][1][k],PI[0][2][k]],[PI[1][0][k],PI[1][1][k],PI[1][2][k]],[PI[2][0][k],PI[2][1][k],PI[2][2][k]]]
        qe=[[Q[0][0][k],Q[0][1][k],Q[0][2][k]],[Q[1][0][k],Q[1][1][k],Q[1][2][k]],[Q[2][0][k],Q[2][1][k],Q[2][2][k]]]
        me=np.dot(np.dot(qe,pe),np.transpose(qe))
        A[k][0]=me[0][0]
        A[k][1]=me[1][1]
        A[k][2]=me[2][2]
        A[k][3]=2*me[0][1]
        A[k][4]=2*me[1][2]
        A[k][5]=2*me[0][2]
        A[k][6]=-2*(me[0][0]*centres[k,0]+me[0][1]*centres[k,1]+me[0][2]*centres[k,2])
        A[k][7]=-2*(me[1][0]*centres[k,0]+me[1][1]*centres[k,1]+me[1][2]*centres[k,2])
        A[k][8]=-2*(me[2][0]*centres[k,0]+me[2][1]*centres[k,1]+me[2][2]*centres[k,2])
        A[k][9]=np.dot(np.dot(centres[k],me),np.transpose(centres[k]))-1

    S=[1,1,1,0,0,0,0,0,0,-RSphere*RSphere]
    return(A,S)

def appendEllipsoidToGroup(sphere, appendFilter, elipspode):
    transform = vtk.vtkTransform()
    print(elipspode)
    transform.Scale(elipspode[0]*2, elipspode[1]*2, elipspode[2]*2)
    transform.Translate(elipspode[3], elipspode[4], elipspode[5])
    transform.RotateX(elipspode[6])
    transform.RotateY(elipspode[7])
    transform.RotateZ(elipspode[8])
    transformer = vtk.vtkTransformPolyDataFilter()
    transformer.SetInputConnection(sphere.GetOutputPort())
    transformer.SetTransform(transform)
    transformer.Update()
    polyData = transformer.GetOutput()


    colors = vtk.vtkUnsignedCharArray()
    colors.SetNumberOfComponents(3)
    colors.SetNumberOfTuples(polyData.GetNumberOfCells())
    r=np.random.randint(256)
    g=np.random.randint(256)
    b=np.random.randint(256)
    for c in range(polyData.GetNumberOfCells()):
        colors.SetTuple(c, [r, g, b])
    polyData.GetCellData().SetScalars(colors)

    if vtk.VTK_MAJOR_VERSION <= 5:
        appendFilter.AddInputConnection(polyData.GetProducerPort())
    else:
        appendFilter.AddInputData(polyData)

def PlotEllipsoideSimple(RSphere, NBEllipsoides, Ellipsoides, pz) :
    sphere = vtk.vtkSphereSource()
    sphere.SetThetaResolution(100)
    sphere.SetPhiResolution(50)
    appendFilter = vtk.vtkAppendPolyData()
    appendEllipsoidToGroup(sphere, appendFilter, [RSphere*2, RSphere*2, RSphere*2, 0, 0, 0, 0, 0, 0])

    for i in range(NBEllipsoides):
        appendEllipsoidToGroup(sphere, appendFilter, Ellipsoides[i])
        appendFilter.Update()

    cleanFilter = vtk.vtkCleanPolyData()
    cleanFilter.SetInputConnection(appendFilter.GetOutputPort())
    cleanFilter.Update()

    # Create a mapper and actor
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(cleanFilter.GetOutputPort())

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
#########################################################################
    plane = vtk.vtkPlane()
    plane.SetOrigin(0, 0, 0)
    plane.SetNormal(0, 1, 0)

    # vtkClipPolyData requires an implicit function to define what it is to
    # clip with. Any implicit function, including complex boolean combinations
    # can be used. Notice that we can specify the value of the implicit function
    # with the SetValue method.
    clipper = vtk.vtkClipPolyData()
    clipper.SetInputConnection(appendFilter.GetOutputPort())
    clipper.SetClipFunction(plane)
    clipper.GenerateClipScalarsOn()
    clipper.GenerateClippedOutputOn()
    clipper.SetValue(0.5)
    clipMapper = vtk.vtkPolyDataMapper()
    clipMapper.SetInputConnection(clipper.GetOutputPort())
    clipMapper.ScalarVisibilityOff()
    backProp = vtk.vtkProperty()
    backProp.SetDiffuseColor(tomato)
    clipActor = vtk.vtkActor()
    clipActor.SetMapper(clipMapper)
    clipActor.GetProperty().SetColor(peacock)
    clipActor.SetBackfaceProperty(backProp)

    # Here we are cutting the cow. Cutting creates lines where the cut
    # function intersects the model. (Clipping removes a portion of the
    # model but the dimension of the data does not change.)
    #
    # The reason we are cutting is to generate a closed polygon at the
    # boundary of the clipping process. The cutter generates line
    # segments, the stripper then puts them together into polylines. We
    # then pull a trick and define polygons using the closed line
    # segments that the stripper created.
    cutEdges = vtk.vtkCutter()
    cutEdges.SetInputConnection(appendFilter.GetOutputPort())
    cutEdges.SetCutFunction(plane)
    cutEdges.GenerateCutScalarsOn()
    cutEdges.SetValue(0, 0.5)
    cutStrips = vtk.vtkStripper()
    cutStrips.SetInputConnection(cutEdges.GetOutputPort())
    cutStrips.Update()
    cutPoly = vtk.vtkPolyData()
    cutPoly.SetPoints(cutStrips.GetOutput().GetPoints())
    cutPoly.SetPolys(cutStrips.GetOutput().GetLines())

    # Triangle filter is robust enough to ignore the duplicate point at
    # the beginning and end of the polygons and triangulate them.
    cutTriangles = vtk.vtkTriangleFilter()
    cutTriangles.SetInputData(cutPoly)
    cutMapper = vtk.vtkPolyDataMapper()
    cutMapper.SetInputData(cutPoly)
    cutMapper.SetInputConnection(cutTriangles.GetOutputPort())
    cutActor = vtk.vtkActor()
    cutActor.SetMapper(cutMapper)

    #cutActor.GetProperty().SetLineWidth(2)
    #cutActor.GetProperty().EdgeVisibilityOn()
    #cutActor.GetProperty().SetColor(peacock)

    # The clipped part of the cow is rendered wireframe.
    restMapper = vtk.vtkPolyDataMapper()
    restMapper.SetInputConnection(clipper.GetClippedOutputPort())
    restMapper.ScalarVisibilityOff()
    restActor = vtk.vtkActor()
    restActor.SetMapper(restMapper)
    restActor.GetProperty().SetRepresentationToWireframe()

    # Create graphics stuff
    ren = vtk.vtkRenderer()
    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(ren)
    iren = vtk.vtkRenderWindowInteractor()

    iren.SetRenderWindow(renWin)

    # Add the actors to the renderer, set the background and size
    ren.AddActor(clipActor)
    ren.AddActor(cutActor)
    ren.AddActor(restActor)
    ren.SetBackground(0, 0, 0)
    ren.ResetCamera()
    ren.GetActiveCamera().Azimuth(30)
    ren.GetActiveCamera().Elevation(30)
    ren.GetActiveCamera().Dolly(1.5)
    ren.ResetCameraClippingRange()

    renWin.SetSize(300, 300)
    iren.Initialize()

    # Lets you move the cut plane back and forth by invoking the function
    # Cut with the appropriate plane value (essentially a distance from
    # the original plane).  This is not used in this code but should give
    # you an idea of how to define a function to do this.
    def Cut(v):
        clipper.SetValue(v)
        cutEdges.SetValue(0, v)
        cutStrips.Update()
        cutPoly.SetPoints(cutStrips.GetOutput().GetPoints())
        cutPoly.SetPolys(cutStrips.GetOutput().GetLines())
        cutMapper.Update()
        renWin.Render()

    def Keypress(obj, event):
        key = obj.GetKeySym()
        if key == "z":
            Cut(1)
        elif key =="a":
            Cut(0)
        elif key =="s":
            Cut(-1)

    iren.AddObserver("KeyPressEvent", Keypress)
    renWin.Render()
    iren.Start()
    return(renWin)

def CloseWin(caller,ev):
    print(caller.GetClassName(), "Event Id:", ev)
    CloseWin.int.GetRenderWindow().Finalize()
    CloseWin.int.TerminateApp()

def GetVolumes(NbEllipsoides,Ellipsoides):
    volumes = 4*np.pi*Ellipsoides[:,0]*Ellipsoides[:,1]*Ellipsoides[:,2]/3
    return(volumes)

def CreateImageFromScene(renWin,filename):
    # Choisir un point de vue, choisir une profondeur de plan et ceer une image corespondant à ce plan
#    extract=vtk.ExtractGrid()
#    extract.SetVOI(xmin,xmax,ymin,ymax,zmin,zmax)
#    extract.SetInput()
#   renWin.GetActiveCamera().Zoom(1.2) #increase zoom
    w2if = vtk.vtkWindowToImageFilter()
    w2if.SetInput(renWin)
    w2if.Update()
    writer = vtk.vtkPNGWriter()
    writer.SetFileName(os.path.join(pathlib.Path(__file__).parent.absolute(),"Images",filename))
    writer.SetInputData(w2if.GetOutput())
    writer.Write()

def ComputeEdgeOfUnion(lnume,A,S,pz):
    use_function_callback = True
    colors = vtk.vtkNamedColors()
    colorsCell=['AliceBlue','Green','Blue','Yellow','Magenta']
    # A renderer and render window
    renderer = vtk.vtkRenderer()
    renderer.SetBackground(colors.GetColor3d('Silver'))

    # render window
    renwin = vtk.vtkRenderWindow()
    renwin.SetWindowName("CallBack")
    renwin.AddRenderer(renderer)

    # An interactor
    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(renwin)

    # Start
    if use_function_callback:
        CloseWin.int=interactor
        interactor.AddObserver('EndInteractionEvent', CloseWin)

    # Create actors
    # create container Sphere
    sphere=vtk.vtkQuadric()
    sphere.SetCoefficients(S)
    sampleS = vtk.vtkSampleFunction()
    sampleS.SetImplicitFunction(sphere)
    # mapper
    mapperSphere = vtk.vtkPolyDataMapper()
    mapperSphere.SetInputConnection(sampleS.GetOutputPort())
    mapperSphere.ScalarVisibilityOff()
#    actorSphere = vtk.vtkActor()
#    actorSphere.SetMapper(mapperSphere)
#    actorSphere.GetProperty().EdgeVisibilityOn()
#    actorSphere.GetProperty().SetColor(colors.GetColor3d('AliceBlue'))
#    actorSphere.GetProperty().SetEdgeColor(colors.GetColor3d('SteelBlue'))
    boundsC=np.array(sampleS.GetModelBounds())
    for k in range(3):
        boundsC[2*k]=-1.1*np.sqrt(-S[9])
        boundsC[2*k+1]=1.1*np.sqrt(-S[9])
    print(S[9])
    print(boundsC)
    sampleS.SetModelBounds(boundsC)

    #    surface.SetValue(0, 0.0)
#
    contoursS = vtk.vtkContourFilter()
    contoursS.SetInputConnection(sampleS.GetOutputPort())
    contoursS.GenerateValues(1, 0.0, 1.2)

    contMapper = vtk.vtkPolyDataMapper()
    contMapper.SetInputConnection(contoursS.GetOutputPort())
    contMapper.SetScalarRange(0.0, 1.2)
    contActor = vtk.vtkActor()
    contActor.GetProperty().SetOpacity(0.1)
    contActor.SetMapper(contMapper)

    # create an ellipsoid using a implicit quadric
    quadric=[None]*len(lnume)
    mapper=[None]*len(lnume)
    actor=[None]*len(lnume)
    sample=[None]*len(lnume)
    contours=[None]*len(lnume)
    bounds=np.array(sampleS.GetModelBounds())
    booleanOperation = vtk.vtkBooleanOperationPolyDataFilter()
    booleanOperation.SetOperationToUnion()
    booleanOperationMapper = vtk.vtkPolyDataMapper()
    booleanOperationActor = vtk.vtkActor()
    for k in lnume:
        quadric[k] = vtk.vtkQuadric()
        quadric[k].SetCoefficients(A[lnume[k]])
        sample[k] = vtk.vtkSampleFunction()
        sample[k].SetImplicitFunction(quadric[k])
        for l in range(3):
            bounds[2*l]=-1.2*np.sqrt(-S[9])
            bounds[2*l+1]=1.2*np.sqrt(-S[9])

        sample[k].SetModelBounds(bounds)
        contours[k] = vtk.vtkContourFilter()
        contours[k].SetInputConnection(sample[k].GetOutputPort())
        # generation d'une valeur dans l'intervalle spécifié (le min et le max sont inclus)
        contours[k].GenerateValues(1, 0.0, 0.1)

        # mapper
        mapper[k] = vtk.vtkPolyDataMapper()
        mapper[k].SetInputConnection(contours[k].GetOutputPort())
        mapper[k].ScalarVisibilityOff()
        actor[k] = vtk.vtkActor()
        actor[k].SetMapper(mapper[k])
        #actor[k].GetProperty().EdgeVisibilityOn()
        actor[k].GetProperty().SetColor(colors.GetColor3d(colorsCell[k%len(colorsCell)]))
        actor[k].GetProperty().SetEdgeColor(colors.GetColor3d('SteelBlue'))
        Ellipsoid2Tri = vtk.vtkTriangleFilter()
        input2 = contours[k].GetOutput()
        Ellipsoid2Tri.SetInputData(input2)
        if k==lnume[0]:
            booleanOperationActor.SetMapper(mapper[lnume[0]])
            EllipsoidsTri= vtk.vtkTriangleFilter()
            input1 = booleanOperation.GetOutputPort().GetOutput() # A revoir
            EllipsoidsTri.SetInputData(input1)
        else :
            booleanOperation.SetInputConnection(0, EllipsoidsTri.GetOutputPort())
            booleanOperation.SetInputConnection(1, Ellipsoid2Tri.GetOutputPort())
            booleanOperation.Update()
            booleanOperationMapper.SetInputConnection(booleanOperation.GetOutputPort())
            booleanOperationMapper.ScalarVisibilityOff()
            booleanOperationActor.GetProperty().SetDiffuseColor(colors.GetColor3d("Banana"))

if __name__ == '__main__':
    main()