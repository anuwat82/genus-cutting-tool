#include "VTK_Header.h"
#include "MouseInteractorStylePP.h"
#include <string>
#include <sstream>

vtkStandardNewMacro(MouseInteractorStylePP);
MouseInteractorStylePP::MouseInteractorStylePP()
{
}
	
void MouseInteractorStylePP::OnLeftButtonDown() 
{
	
	

		vtkPointPicker* picker = static_cast<vtkPointPicker *>(this->Interactor->GetPicker());
		//picker->InitializePickList(); 	
		
		int ret = picker->Pick(	this->Interactor->GetEventPosition()[0], 
						this->Interactor->GetEventPosition()[1], 
						0,  // always zero.
						this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer());
		
		
	
		vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
	
}
