#pragma once
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkObjectFactory.h>
#include <vtkCommand.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderWindow.h>
#include <vtkRendererCollection.h>
#include <vtkPointPicker.h>
#include <vtkCellPicker.h>
class MouseInteractorStylePP :  public vtkInteractorStyleTrackballCamera
{
public:
	vtkTypeMacro(MouseInteractorStylePP, vtkInteractorStyleTrackballCamera);
	static MouseInteractorStylePP* New();
    

	virtual void OnLeftButtonDown();
	
protected: 
	MouseInteractorStylePP();
protected:
	
};


 
    

