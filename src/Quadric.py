import sys
import vtk
import numpy as np
import math
import os
import pathlib
import subprocess

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

                (A,S)=FromGeomToCartForEllipsoid(NbEllipsoides,Ellipsoides,RSphere)
                lnume=np.array(range(NbEllipsoides))
                volumes=GetVolumes(NbEllipsoides,Ellipsoides)
                #print(volumes)

                # Affichage es ellipsoides et de la sphère englobante
                renWin=PlotEllipsoideSimple(lnume,A,S,pz)

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

def PlotEllipsoideSimple(lnume,A,S,pz) :
    # lnume : liste des ellipsoides
    # A : coefficients des equations pour les ellipsoides
    # S : coefficients des eéquations pour la sphère anglobante
    # pz : pourcentage entre -& et 1 pour positionner le plan de coupe a une altitude pz*Rayon


    use_function_callback = False

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
    bounds=np.array(sampleS.GetModelBounds())
    for k in range(3):
        bounds[2*k]=-1.1*np.sqrt(-S[9])
        bounds[2*k+1]=1.1*np.sqrt(-S[9])
    print(S[9])
    print(bounds)
    sampleS.SetModelBounds(bounds)
#    # contour
#    surface = vtk.vtkContourFilter()
#    surface.SetInputConnection(sample.GetOutputPort())
##    surface.SetInputConnection(extract.GetOutputPort())
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
     #Add the actors to the scene
    renderer.AddActor(contActor)

    for k in lnume:
        renderer.AddActor(actor[k])

    interactor.Initialize()
    renwin.Render()
    interactor.Start()
    return(renwin)

    colors = vtk.vtkNamedColors()

    #Set the background color.
    colors.SetColor("BkgColor", [26, 51, 77, 255])

    # Create a plane
    planeSource = vtk.vtkPlaneSource()
    planeSource.SetCenter(1.0, 0.0, 0.0)
    planeSource.SetNormal(1.0, 0.0, 1.0)
    planeSource.Update()

    plane = planeSource.GetOutput()

    # Create a mapper and actor
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(plane)

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(colors.GetColor3d("Cyan"))

    # Create a renderer, render window and interactor
    renderer = vtk.vtkRenderer()
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.SetWindowName("Plane")
    renderWindow.AddRenderer(renderer)
    renderWindowInteractor = vtk.vtkRenderWindowInteractor()
    renderWindowInteractor.SetRenderWindow(renderWindow)

    # Add the actors to the scene
    renderer.AddActor(actor)
    renderer.SetBackground(colors.GetColor3d("BkgColor"))

    # Render and interact
    renderWindow.Render()
    renderWindowInteractor.Start()

def PlotEllipsoide(lnume,A,S,pz) :
    # lnume : liste des ellipsoides
    # A : coefficients des equations pour les ellipsoides
    # S : coefficients des eéquations pour la sphère anglobante
    # pz : pourcentage entre -& et 1 pour positionner le plan de coupe a une altitude pz*Rayon


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
    bounds=np.array(sampleS.GetModelBounds())
    for k in range(3):
        bounds[2*k]=-1.1*np.sqrt(-S[9])
        bounds[2*k+1]=1.1*np.sqrt(-S[9])
    print(S[9])
    print(bounds)
    sampleS.SetModelBounds(bounds)
#    # contour
#    surface = vtk.vtkContourFilter()
#    surface.SetInputConnection(sample.GetOutputPort())
##    surface.SetInputConnection(extract.GetOutputPort())
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


    #Add the actors to the scene
    renderer.AddActor(contActor)

    for k in lnume:
        renderer.AddActor(actor[k])

## Add cut of sphere
#
#    cutterS = vtk.vtkCutter()
#    cutterS.SetInputConnection(sampleS.GetOutputPort())
#    plane = vtk.vtkPlane()
#
#    #Plan de coupe xy
#    altitude=pz*np.sqrt(-S[9])
#    plane.SetOrigin([0,0,altitude])
#    plane.SetNormal([0,0,1])
#
#    cutterS.SetCutFunction(plane)
##    cutterS.SetInputConnection(sampleS.GetOutputPort())
##    cutterS.GenerateCutScalarsOn()
##    cutterS.SetValue(0,0.5)
##    cutterS.Update()
#
#    cutStrips = vtk.vtkStripper()
#    cutStrips.SetInputConnection(cutterS.GetOutputPort())
#    cutStrips.Update()
#
#    circle= cutStrips.GetOutput()
#
#    # write circle out
#    polyDataWriter = vtk.vtkXMLPolyDataWriter()
#    if vtk.VTK_MAJOR_VERSION <= 5:
#        polyDataWriter.SetInput(circle)
#    else:
#        polyDataWriter.SetInputData(circle)
#
#    polyDataWriter.SetFileName("circle.vtp")
#    polyDataWriter.SetCompressorTypeToNone()
#    polyDataWriter.SetDataModeToAscii()
#    polyDataWriter.Write()
#
## prepare the binary image's voxel grid
#    whiteImage = vtk.vtkImageData()
#    boundsI = [0]*6
#    circle.GetBounds(boundsI)
#    spacing = [0]*3 # desired volume spacing
#    spacing[0] = 0.5
#    spacing[1] = 0.5
#    spacing[2] = 0.5
#    whiteImage.SetSpacing(spacing)
#
## compute dimensions
#    dim = [0]*3
#    for i in range(3):
#        dim[i] = int(math.ceil((bounds[i * 2 + 1] - bounds[i * 2]) / spacing[i])) + 1
#        if (dim[i] < 1):
#            dim[i] = 1
#    whiteImage.SetDimensions(dim)
#    whiteImage.SetExtent(0, dim[0] - 1, 0, dim[1] - 1, 0, dim[2] - 1)
#    origin = [0]*3
## NOTE: I am not sure whether or not we had to add some offset!
#    origin[0] = bounds[0]# + spacing[0] / 2
#    origin[1] = bounds[2]# + spacing[1] / 2
#    origin[2] = bounds[4]# + spacing[2] / 2
#    whiteImage.SetOrigin(origin)
#    if vtk.VTK_MAJOR_VERSION <= 5:
#        whiteImage.SetScalarTypeToUnsignedChar()
#        whiteImage.AllocateScalars()
#    else:
#        whiteImage.AllocateScalars(vtk.VTK_UNSIGNED_CHAR, 1)
#
## fill the image with foreground voxels:
#    inval = 255
#    outval = 0
#    count = whiteImage.GetNumberOfPoints()
##for (vtkIdType i = 0 i < count ++i)
#    for i in range(count):
#        whiteImage.GetPointData().GetScalars().SetTuple1(i, inval)
#
## sweep polygonal data (this is the important thing with contours!)
#    extruder = vtk.vtkLinearExtrusionFilter()
#    if vtk.VTK_MAJOR_VERSION <= 5:
#        extruder.SetInput(circle)
#    else:
#        extruder.SetInputData(circle)
#
#    extruder.SetScaleFactor(1.)
#    extruder.SetExtrusionTypeToNormalExtrusion()
#    extruder.SetVector(0, 0, 1)
#    extruder.Update()
#
#    # polygonal data -. image stencil:
#    pol2stenc = vtk.vtkPolyDataToImageStencil()
#    pol2stenc.SetTolerance(0) # important if extruder.SetVector(0, 0, 1) !!!
#    pol2stenc.SetInputConnection(extruder.GetOutputPort())
#    pol2stenc.SetOutputOrigin(origin)
#    pol2stenc.SetOutputSpacing(spacing)
#    pol2stenc.SetOutputWholeExtent(whiteImage.GetExtent())
#    pol2stenc.Update()
#
## cut the corresponding white image and set the background:
#    imgstenc = vtk.vtkImageStencil()
#    if vtk.VTK_MAJOR_VERSION <= 5:
#        imgstenc.SetInput(whiteImage)
#        imgstenc.SetStencil(pol2stenc.GetOutput())
#    else:
#        imgstenc.SetInputData(whiteImage)
#        imgstenc.SetStencilConnection(pol2stenc.GetOutputPort())
#
#    imgstenc.ReverseStencilOff()
#    imgstenc.SetBackgroundValue(outval)
#    imgstenc.Update()
#
#    imageWriter = vtk.vtkMetaImageWriter()
#    imageWriter.SetFileName("labelImage.mhd")
#    imageWriter.SetInputConnection(imgstenc.GetOutputPort())
#    imageWriter.Write()
#
#    imageWriter = vtk.vtkPNGWriter()
#    imageWriter.SetFileName("labelImage.png")
#    imageWriter.SetInputConnection(imgstenc.GetOutputPort())
#    imageWriter.Write()
#
#
#    cutPoly = vtk.vtkPolyData()
#
#    cutPoly.SetPoints(cutStrips.GetOutput().GetPoints())
#    cutPoly.SetPolys(cutStrips.GetOutput().GetLines())
#
#    # Get the edges from the mesh
#    edges = vtk.vtkExtractEdges()
#    edges.SetInput(cutPoly)
#    edge_mapper = vtk.vtkPolyDataMapper()
#    edge_mapper.SetInput(edges.GetOutput())
#
## Make an actor for those edges
#    edge_actor = vtk.vtkActor()
#    edge_actor.SetMapper(edge_mapper)
## Make the actor red (there are other ways of doing this also)
#    edge_actor.GetProperty().SetColor(1,0,0)
#
#    renderer.AddActor(edge_actor)
#
##    # Avoid z-buffer fighting
#    vtk.vtkPolyDataMapper().SetResolveCoincidentTopologyToPolygonOffset()
##
##    # Extract boundary from cutPoly
##    cutBoundary = vtk.vtkFeatureEdges()
##    cutBoundary.SetInputConnection(cutPoly)
##
##    cutBoundary.Update()
#    cutterSMapper = vtk.vtkPolyDataMapper()
#    cutterSMapper.SetInputConnection(cutterS.GetOutputPort())
#    cutterSMapper.ScalarVisibilityOff()
##    cutterSMapper.SetScalarRange(0.,0.1)
###    cutterMapper.SetLookupTable(clut)
###
#    cut= vtk.vtkActor()
#    cut.SetMapper(cutterSMapper)
#    cut.GetProperty().SetOpacity(1)
#    cut.GetProperty().EdgeVisibilityOn()
#    cut.GetProperty().SetColor( 1,0,0 )
#    cut.GetProperty().SetEdgeColor( 1,1,1 )
#
#    renderer.AddActor(cut)
#
##    ren1.AddViewProp( triActor )
##    ren1.SetViewport( 0,0,0.5,1.0)
#
#    cutters=[None]*len(lnume)
#    cuttersMapper=[None]*len(lnume)
#    cut=[None]*len(lnume)
#    for k in lnume:
#         cutters[k] = vtk.vtkCutter()
#         cutters[k].SetCutFunction(plane)
#         cutters[k].SetInputConnection(sample[k].GetOutputPort())
#         cutters[k].Update()
#         cuttersMapper[k] = vtk.vtkPolyDataMapper()
#         cuttersMapper[k].SetInputConnection(cutters[k].GetOutputPort())
#         cuttersMapper[k].SetScalarRange(0.,0.1)
#         #    cutterMapper.SetLookupTable(clut)
#         #
#         cut[k]= vtk.vtkActor()
#         cut[k].SetMapper(cuttersMapper[k])
#         cut[k].GetProperty().SetOpacity(1)
##         renderer.AddActor(cut[k])
#
#    # The sample function generates a distance function from the implicit
#    # function. This is then contoured to get a polygonal surface.
##    sample = vtk.vtkSampleFunction()
##    sample.SetImplicitFunction(quadric)
##
##    booleanUnion = vtk.vtkImplicitBoolean()
##    booleanUnion.AddFunction(sphere)
##    booleanUnion.SetOperationType(0)
##    extract = vtk.vtkExtractGeometry()
##    extract.SetInputConnection(sample.GetOutputPort())
##    extract.SetImplicitFunction(booleanUnion)
##    #sample.SetImplicitFunction(sphere)
##    #print(sample.GetModelBounds())
##    #sample.SetModelBounds(-0.5, 0.5, -0.5, 0.5, -0.5, 0.5)
##    bounds=np.array(sample.GetModelBounds())
##    for k in range(3):
##        bounds[2*k]=10*np.sqrt(1/A[nume][k])*bounds[2*k]
##        bounds[2*k+1]=10*np.sqrt(1/A[nume][k])*bounds[2*k+1]
##
##    #for k in range(3):
##     #   bounds[2*k]=np.sqrt(-S[9]/S[k])*bounds[2*k]
##      #  bounds[2*k+1]=np.sqrt(-S[9]/S[k])*bounds[2*k+1]
##
##    sample.SetModelBounds(bounds)
##    print(sample.GetModelBounds())
##
##
##    sample.SetSampleDimensions(40, 40, 40)
##    sample.ComputeNormalsOff()
##
##    # contour
##    surface = vtk.vtkContourFilter()
##    #surface.SetInputConnection(sample.GetOutputPort())
##    surface.SetInputConnection(extract.GetOutputPort())
##    surface.SetValue(0, 0.0)
##
##
##
##
##
##    # add the actor
##    renderer.AddActor(actor)
#
#
#

    interactor.Initialize()
    renwin.Render()
    interactor.Start()
    return(renwin)

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
    writer
    writer.SetInputData(w2if.GetOutput())
    writer.Write()

    ### Translate and Rotate

#    center_x, center_y, center_z = actor.GetCenter()
#    ren.AddActor(actor)
#    w = vtk.vtkTransform()
#    w.Translate(-center_x, -center_y, -center_z)
#    w.RotateX(coordinates[0])
#    w.RotateZ(coordinates[2])
#    actor.SetUserTransform(w)


#    plane=vtk.vtkPlane()
#    plane.SetOrigin(0,4,2)
#    plane.SetNormal(0,0,1)
#

#
#    cutter=vtk.vtkCutter()
#    cutter.SetInput(extract.GetOutput())
#    cutter.SetCutFunction(plane)
#    cutter.GenerateCutScalarsOff()
#    cutter.SetSortByToSortCellByCell()
#
#    clut= vtk.vtkLookupTable()
#    clut.SetHueRange(0,0.67)
#    clut.Build()
#
#    cutterMapper = vtk.vtkPolyDataMapper()
#    cutterMapper.SetInput(cutter.GetOutput())
#    cutterMapper.SetScalarRange(.18,.7)
#    cutterMapper.SetLookupTable(clut)
#
#    cut= vtk.vtkActor()
#    cut.SetMapper(cutterMapper)
#    cut.GetProperty().SetOpacity(1)
#
#    cast=vtkImageCast()
#    cast.SetInput()
#    cast.SetOutputScalartTypeToUnsignedChar()

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

# intersection entre scene3d et les ellip
