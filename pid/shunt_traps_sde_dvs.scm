; *** Al-BSF solar cell for Na migration***
; Guillaume von Gastrow
; These command are written in Scheme programming language (based on LISP)

; *** INITIALIZATION ***
; clear structure
(sde:clear)

; MATERIALS
(define CellMaterial 			"Silicon")
(define ARCMaterial 			"Si3N4")
(define shuntmat 					"Metal")
(define contactWidth 			${contact_length})
(define cellWidth 				${device_length})
(define THICKNESS 				200)
(define emitterThickness 	0.6)
; The thickness of the Al BSF
(define tBSF 			0.3)
; An error flag TODO: Consider removing it
(define err_flag 	0)
; Silicon nitride thickness
(define tSiNx 		75e-3)
; Shunt width
(define shuntw 		0.01) ; um
; Shunt depths
(define ${dshunt_name} ${shunt_depth}) ; modified by the batch Sentaurus script based on the Na profile

; New-replace-old option (default)
(sdegeo:set-default-boolean "ABA") ; subtract overlapping regions from the existing regions

; *** CHOICE OF DOPING PROFILE ***
;(define dopingmodel "uniform")
(define dopingmodel "gauss")

; DOPING PARAMETERS
(define em_doping 1e19) ; Uniform doping value or surface dopant concentration, depending on doping model
(define base_doping 1e16) ; Silicon wafer base doping

; MESH PARAMETERS
(define xmax 10)
(define xmin 0.1)
(define ymax 10)
(define ymin 0.1)
; Import optical generation file
(define OptGenFile "${optical_generation_file}" )
; Define the number of shunts
(define shuntBoxStart (list ${shunt_positions}) )
(define shuntBoxNumber (list ${shunt_index}) )
(define shuntDY ${shunt_dy})

; *** GEOMETRY ***
; convention: x=length y=thickness

; create SiNx film (leave space for the metal contact before)
(sdegeo:create-rectangle
	(position (/ contactWidth 2) 0 0)
	(position (- cellWidth (/ contactWidth 2)) (- tSiNx) 0)
	ARCMaterial
	"SiNx-region"
)

; Choose between uniform doping and error function doping (ifelse statement, similar to switch syntax)
(cond
  ((string=? dopingmodel "uniform") ; condition 1, uniform doping
		(begin
			(newline)
			(display "Using an uniform doping profile")
			(newline)
			; create n-emitter
			(sdegeo:create-rectangle (position 0 0 0) (position cellWidth emitterThickness 0) CellMaterial "Si-profile-region")
			; create p-emitter base
			(sdegeo:create-rectangle (position 0 emitterThickness 0) (position cellWidth THICKNESS 0) CellMaterial "Si-base-region")

			; *** DOPING ****
			; Emitter
			(sdedr:define-constant-profile "emitter-profile" "PhosphorusActiveConcentration" em_doping)
			(sdedr:define-constant-profile-region "em-doping-placement" "emitter-profile" "Si-profile-region")

			; Base
			(sdedr:define-constant-profile "base-doping-profile" "BoronActiveConcentration" base_doping)
			(sdedr:define-constant-profile-region "base-doping-placement" "base-doping-profile" "Si-base-region") ; place the base doping profile in the base region
		) ; end of begin
  ) ; end of ((string=? dopingmodel "uniform")
  ((string=? dopingmodel "gauss")
	(begin
		(display "Using an error function doping profile")
		(newline)
		; create Si region
		(sdegeo:create-rectangle (position 0 0 0) (position cellWidth THICKNESS 0) CellMaterial "Si-profile-region") ; unique region on which an analytical doping profile will be added
		; Base
		(sdedr:define-constant-profile "base-doping-profile" "BoronActiveConcentration" base_doping)
		(sdedr:define-constant-profile-region "base-doping-placement" "base-doping-profile" "Si-profile-region") ; place the base doping profile in the base region
		; *** DOPING ****
		; Emitter
		(sdedr:define-gaussian-profile "Gauss.Emitter"
			"PhosphorusActiveConcentration"
			"PeakPos" 0.0
			"PeakVal" em_doping
			"ValueAtDepth" base_doping
			"Depth" emitterThickness
			"Gauss" "Factor" 0.5) ; error function doping profile. C0 is "MaxVal" (surface conc is C0/2) and "Length" is 2sqrt(Dt). See sde manual p. 222 and smesh manual p.114. Parameters based on profiles from Kerr et al in Fig. 4, JAP, 89, 2001.
		(sdedr:define-refinement-window "RefEvalWin.Emitter"
			"Line"
			(position 0 0 0)
			(position cellWidth 0 0)
		) ; window for the analytical profile placement. In 2D, must be a line normal to the profile.
		(sdedr:define-analytical-profile-placement "PlaceAP.Emitter"
			"Gauss.Emitter"
			"RefEvalWin.Emitter"
			"Positive"
			"NoReplace"
			"Eval"
		) ; Place the analytical profile-doping-window") ; not used because base-doping-window not defined
	)
  )
) ; end of cond


(sdedr:define-gaussian-profile "Gauss.AlBSF"
	"AluminumActiveConcentration"
	"PeakPos" 0
	"PeakVal" em_doping
	"ValueAtDepth" base_doping
	"Depth" 1
	"Gauss" "Factor" 0.8)
(sdedr:define-refinement-window "RefEvalWin.AlBSF"
	"Line"
	(position 0 THICKNESS 0)
	(position cellWidth THICKNESS 0)
) ; window for the analytical profile placement.
(sdedr:define-analytical-profile-placement "PlaceAP.AlBSF"
	"Gauss.AlBSF"
	"RefEvalWin.AlBSF"
	"Negative"
	"NoReplace" "Eval"
) ; Place the analytical profile


; *** CONTACTS ***
(sdegeo:define-contact-set "em_contact" 	4 (color:rgb 1 0 0) "##") 	; em_contact
(sdegeo:define-contact-set "base_contact" 4 (color:rgb 0 0 1) "##") 	; base_contact
; Create the metal region
(sdegeo:create-rectangle (position 0 0 0) (position (/ contactWidth 2) (- 0.3) 0) "Silver" "EmitterContact")
(sdegeo:create-rectangle (position (- cellWidth (/ contactWidth 2)) 0 0) (position cellWidth (- 0.3) 0) "Silver" "EmitterContact")
(sdegeo:delete-region (find-body-id (position (/ contactWidth 4) (- 0.3) 0)))
(sdegeo:delete-region (find-body-id (position (- cellWidth (/ contactWidth 4)) (- 0.3) 0)))

(sdegeo:create-rectangle (position 0 THICKNESS 0) (position cellWidth (+ THICKNESS 0.3) 0) "Aluminum" "Region.BSF")
(sdegeo:delete-region (find-body-id (position (/ cellWidth 2) (+ THICKNESS 0.15) 0) ))

; Activate contacts
(sdegeo:set-current-contact-set "em_contact")
(sdegeo:define-2d-contact (find-edge-id (position (/ contactWidth 4) 0 0)) "em_contact")
(sdegeo:define-2d-contact (find-edge-id (position (- cellWidth (/ contactWidth 4)) 0 0)) "em_contact")

(sdegeo:set-current-contact-set "base_contact")
(sdegeo:define-2d-contact (find-edge-id (position (/ cellWidth 2) THICKNESS 0)) "base_contact")



; *** MESH ***
; * WHOLE DOMAIN
(sdedr:define-refeval-window "domain-ref" "Rectangle" (position 0 (- tSiNx) 0) (position cellWidth THICKNESS 0))
(sdedr:define-refinement-size "domain-ref-size" xmax ymax xmin ymin)
(sdedr:define-refinement-placement "domain-ref-pl" "domain-ref-size" "domain-ref")

(sdedr:define-refeval-window "RefWin.All"
	"Rectangle"
	(position 0 0 0)
	(position cellWidth THICKNESS 0)
)
(sdedr:define-refinement-size "RefDef.All"
	(/ xmax 20) (/ ymax 20)
	(/ xmin 25) (/ ymin 25)
)
(sdedr:define-refinement-function "RefDef.All"
	"DopingConcentration"
	"MaxTransDiff" 1
)

(sdedr:define-refinement-function "RefDef.All"
	"PhosphorusActiveConcentration"
	"MaxTransDiff" 1
)

(sdedr:define-refinement-function "RefDef.All"
	"AluminumActiveConcentration"
	"MaxTransDiff" 1
)

(sdedr:define-refinement-placement "PlaceRF.all" "RefDef.All" "RefWin.All")



; * p-n JUNCTION REFINEMENT * Used for an abrupt junction but shouldn't be necessary for a erf doping profile
(if (string=? dopingmodel "uniform")
	(begin
		(sdedr:define-refeval-window "junction-ref" "Rectangle" (position 0 (- emitterThickness 0.050) 0) (position cellWidth (+ emitterThickness 0.050) 0))
		(sdedr:define-refinement-size "junction-ref-size" (/ xmax 10) (/ ymax 10) (/ xmin 10) (/ ymin 10))
		(sdedr:define-refinement-placement "junction-ref-pl" "junction-ref-size" "junction-ref")
	)
)

; Mesh refinement in the doped region (go up to 1 um)
(sdedr:define-refeval-window "junction-ref" "Rectangle" (position 0 0 0) (position cellWidth 1.2 0))
(sdedr:define-refinement-size "junction-ref-size" (/ xmax 20) (/ ymax 20) (/ xmin 50) (/ ymin 50))
(sdedr:define-refinement-placement "junction-ref-pl" "junction-ref-size" "junction-ref")



(display "  front contact refinement") (newline)
(sdedr:define-refinement-window "RefWin.FrontContact1" "Rectangle"
	(position  0 0 0)
	(position  (+ (/ contactWidth 2) 0.5) THICKNESS 0)
)
(sdedr:define-refinement-window "RefWin.FrontContact2" "Rectangle"
	(position  (- cellWidth (- (/ contactWidth 2 ) 0.5)) 0 0)
	(position  cellWidth THICKNESS 0)
)

(sdedr:define-refinement-size "RefDef.FrontContact"
	(/ contactWidth 10) (/ THICKNESS 10)
	(/ contactWidth 20) (/ THICKNESS 50)
)

(sdedr:define-refinement-placement "PlaceRF.FrontContact1" "RefDef.FrontContact" "RefWin.FrontContact1" )
(sdedr:define-refinement-placement "PlaceRF.FrontContact2" "RefDef.FrontContact" "RefWin.FrontContact2" )

; *** OPTICAL GENERATION PROFILE ***
; NOTE: Make sure the window is defined in a direction normal to the optical generation profile!
(sdedr:define-refinement-window "opt_win" "Line" (position (* contactWidth 0.5) 0 0) (position (- cellWidth (* contactWidth 0.5)) 0 0))
(sdedr:define-1d-external-profile "1d_opt_def" OptGenFile "Scale" 1.0 "Range" 0 THICKNESS "Erf" "Factor" 0) ;
(sdedr:define-analytical-profile-placement "opt_place" "1d_opt_def" "opt_win" "Positive" "NoReplace" "Eval")


; * BUILDING MESH
(sde:build-mesh "snmesh" "-a -c boxmethod" "${nodenum}") ; save mesh files with name "nodnum"
