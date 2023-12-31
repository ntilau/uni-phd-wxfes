/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkProp3DButtonRepresentation.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
/**
 * @class   vtkProp3DButtonRepresentation
 * @brief   defines a representation for a vtkButtonWidget
 *
 * This class implements one type of vtkButtonRepresentation. Each button
 * state can be represented with a separate instance of vtkProp3D. Thus
 * buttons can be represented with vtkActor, vtkImageActor, volumes (e.g.,
 * vtkVolume) and/or any other vtkProp3D. Also, the class invokes events when
 * highlighting occurs (i.e., hovering, selecting) so that appropriate action
 * can be taken to highlight the button (if desired).
 *
 * To use this representation, always begin by specifying the number of
 * button states.  Then provide, for each state, an instance of vtkProp3D.
 *
 * This widget representation uses the conventional placement method. The
 * button is placed inside the bounding box defined by PlaceWidget by translating
 * and scaling the vtkProp3D to fit (each vtkProp3D is transformed). Therefore,
 * you must define the number of button states and each state (i.e., vtkProp3D)
 * prior to calling vtkPlaceWidget.
 *
 * @sa
 * vtkButtonWidget vtkButtonRepresentation vtkButtonSource vtkEllipticalButtonSource
 * vtkRectangularButtonSource
 */

#ifndef vtkProp3DButtonRepresentation_h
#define vtkProp3DButtonRepresentation_h

#include "vtkButtonRepresentation.h"
#include "vtkInteractionWidgetsModule.h" // For export macro

class vtkPropPicker;
class vtkProp3D;
class vtkProp3DFollower;
class vtkPropArray; // PIMPLd

class VTKINTERACTIONWIDGETS_EXPORT vtkProp3DButtonRepresentation : public vtkButtonRepresentation
{
public:
  /**
   * Instantiate the class.
   */
  static vtkProp3DButtonRepresentation* New();

  //@{
  /**
   * Standard methods for instances of the class.
   */
  vtkTypeMacro(vtkProp3DButtonRepresentation, vtkButtonRepresentation);
  void PrintSelf(ostream& os, vtkIndent indent) override;
  //@}

  //@{
  /**
   * Add the ith texture corresponding to the ith button state.
   * The parameter i should be (0 <= i < NumberOfStates).
   */
  void SetButtonProp(int i, vtkProp3D* prop);
  vtkProp3D* GetButtonProp(int i);
  //@}

  //@{
  /**
   * Specify whether the button should always face the camera. If enabled,
   * the button reorients itself towards the camera as the camera moves.
   */
  vtkSetMacro(FollowCamera, vtkTypeBool);
  vtkGetMacro(FollowCamera, vtkTypeBool);
  vtkBooleanMacro(FollowCamera, vtkTypeBool);
  //@}

  /**
   * Extend the vtkButtonRepresentation::SetState() method.
   */
  void SetState(int state) override;

  //@{
  /**
   * Provide the necessary methods to satisfy the vtkWidgetRepresentation API.
   */
  int ComputeInteractionState(int X, int Y, int modify = 0) override;
  void BuildRepresentation() override;
  //@}

  /**
   * This method positions (translates and scales the props) into the
   * bounding box specified. Note all the button props are scaled.
   */
  void PlaceWidget(double bounds[6]) override;

  //@{
  /**
   * Provide the necessary methods to satisfy the rendering API.
   */
  void ShallowCopy(vtkProp* prop) override;
  double* GetBounds() override;
  void GetActors(vtkPropCollection* pc) override;
  void ReleaseGraphicsResources(vtkWindow*) override;
  int RenderOpaqueGeometry(vtkViewport*) override;
  int RenderVolumetricGeometry(vtkViewport*) override;
  int RenderTranslucentPolygonalGeometry(vtkViewport*) override;
  vtkTypeBool HasTranslucentPolygonalGeometry() override;
  //@}

  /*
   * Register internal Pickers within PickingManager
   */
  void RegisterPickers() override;

protected:
  vtkProp3DButtonRepresentation();
  ~vtkProp3DButtonRepresentation() override;

  // The current vtkProp3D used to represent the button
  vtkProp3D* CurrentProp;

  // Follow the camera if requested
  vtkProp3DFollower* Follower;
  vtkTypeBool FollowCamera;

  // Keep track of the props associated with the N
  // states of the button. This is a PIMPLd stl map.
  vtkPropArray* PropArray;

  // For picking the button
  vtkPropPicker* Picker;

private:
  vtkProp3DButtonRepresentation(const vtkProp3DButtonRepresentation&) = delete;
  void operator=(const vtkProp3DButtonRepresentation&) = delete;
};

#endif
