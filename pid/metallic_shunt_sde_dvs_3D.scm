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
(define cellLength 				25)
(define THICKNESS 				200)
(define emitterThickness 	0.6)
(define thicknessARC 			75e-3)
(define shuntWidth 				0.01) ; um
(define shuntLength				0.01) ; um
; The thickness of the Al BSF
(define tBSF 			0.3)
; An error flag TODO: Consider removing it
(define err_flag 	0)
; Shunt depths
(define ${dshunt_name} ${shunt_depth}) ; modified by the batch Sentaurus script based on the Na profile

; New-replace-old option (default)
(sdegeo:set-default-boolean "ABA") ; subtract overlapping regions from the existing regions

; *** CHOICE OF DOPING PROFILE ***
;(define dopingmodel "uniform")
(define dopingmodel "gauss")

; *** DEFINITIONS ***
; Initialize error file errfile.txt to no error (0)
(call-with-output-file "errfile.txt"
	(lambda (file)
		(write 0 file)))

; Position of shunt1
(define shunt1_pos_x1 (+ contactWidth (/ (- cellWidth contactWidth) 2)) ) ; left limit on the x-axis
(define shunt1_pos_x2 (+ (+ contactWidth (/ (- cellWidth contactWidth) 2) ) shuntWidth) ) ; right limit on the x-axis
(define shunt1_pos_y1 (+ contactWidth (/ cellLength 2) ) )
(define shunt1_pos_y2 (+ (+ contactWidth (/ cellLength 2) ) shuntLength) )

(display shunt1_pos_x1)

; DOPING PARAMETERS
(define em_doping 1e19) ; Uniform doping value or surface dopant concentration, depending on doping model
(define base_doping 1e16) ; Silicon wafer base doping

; MESH PARAMETERS
(define xmax 10)
(define xmin 0.1)
(define ymax 5)
(define ymin 0.1)
(define zmax 10)
(define zmin 0.1)
; Import optical generation file
(define OptGenFile "${optical_generation_file}" )
; Define the number of shunts
(define shuntBoxStart (list ${shunt_positions}) )
; Define the index of the shunts to construct the names that sdevice refers to in the region-wise definition of Resist0
(define shuntBoxNumber (list ${shunt_index}) )
(define shuntDY ${shunt_dy})

; *** GEOMETRY ***
; convention: x=length y=thickness
; create SiNx film (leave space for the metal contact before)
(sdegeo:create-cuboid (position contactWidth 0 0) (position cellWidth cellLength thicknessARC) ARCMaterial "SiNx-region")

; Choose between uniform doping and error function doping (ifelse statement, similar to switch syntax)
(cond
  ((string=? dopingmodel "uniform") ; condition 1, uniform doping
		(begin
			(newline)
			(display "Using an uniform doping profile")
			(newline)
			; create n-emitter
			(sdegeo:create-cuboid (position 0 0 0) (position cellWidth  cellLength (- emitterThickness) ) CellMaterial "Si-profile-region")
			; create p-emitter base
			(sdegeo:create-cuboid (position 0 emitterThickness 0) (position cellWidth cellLength (- THICKNESS) ) CellMaterial "Si-base-region")

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
		(sdegeo:create-cuboid (position 0 0 0) (position cellWidth cellLength (- THICKNESS) ) CellMaterial "Si-profile-region") ; unique region on which an analytical doping profile will be added
		; Base
		(sdedr:define-constant-profile "base-doping-profile" "BoronActiveConcentration" base_doping)
		(sdedr:define-constant-profile-region "base-doping-placement" "base-doping-profile" "Si-profile-region") ; place the base doping profile in the base region

		; *** DOPING ****
		; Emitter
		;(sdedr:define-constant-profile "emitter-profile" "PhosphorusActiveConcentration" em_doping) ; constant emitter doping
		; (sdedr:define-erf-profile "emitter-profile" "BoronActiveConcentration" "SymPos" 0.0 "MaxVal" (* 2 em_doping) "Length" 0.2 "erf" "factor" 1) ; error function doping profile. C0 is "MaxVal" (surface conc is C0/2) and "Length" is 2sqrt(Dt). See sde manual p. 222 and smesh manual p.114. Parameters based on profiles from Kerr et al in Fig. 4, JAP, 89, 2001.
		(sdedr:define-gaussian-profile "Gauss.Emitter"
			"PhosphorusActiveConcentration"
			"PeakPos" 0.0
			"PeakVal" em_doping
			"ValueAtDepth" base_doping
			"Depth" emitterThickness
			"Gauss" "Factor" 0.8) ; error function doping profile. C0 is "MaxVal" (surface conc is C0/2) and "Length" is 2sqrt(Dt). See sde manual p. 222 and smesh manual p.114. Parameters based on profiles from Kerr et al in Fig. 4, JAP, 89, 2001.
		(sdedr:define-refinement-window "RefEvalWin.Emitter"
			"Rectangle"
			(position 0 0 0)
			(position cellWidth cellLength 0)
		) ; window for the analytical profile placement. In 2D, must be a line normal to the profile.
		(sdedr:define-analytical-profile-placement "PlaceAP.Emitter"
			"Gauss.Emitter"
			"RefEvalWin.Emitter"
			"Negative"
			"NoReplace"
			"Eval"
		) ; Place the analytical profile
	)
  )
) ; end of cond

; *** SHUNTS ***
; create shunts, ie conductive regions through the PN junction with recombination centers at the interface with silicon
; Later will need to place shunts with a do-loop across the width of the PN junction
; position of the shunt contactWidth+(cellWidth-contactWidth)/2

(display ${dshunt_name})
(display "\n")

; Define the shunts
(if (> shuntDY 0)
	(begin
		(for-each
			(lambda (SHUNTY INDEX)
				(begin
					(define REGION (string-append "shunt.Region." (number->string INDEX) ))
					(sdegeo:create-cuboid
						(position shunt1_pos_x1 shunt1_pos_y1 (- SHUNTY) )
						(position shunt1_pos_x2 shunt1_pos_y2 (- (+ SHUNTY shuntDY)) )
						shuntmat
						REGION
					) ; end of sdegeo:create-rectangle
					(display "Created region ")
					(display REGION)
					(newline)
				) ; end of begin
			) ; end of lambda
			shuntBoxStart shuntBoxNumber
		) ; end of for-each
	) ; end of begin
) ; end of if

(sdedr:define-gaussian-profile "Gauss.AlBSF"
	"AluminumActiveConcentration"
	"PeakPos" 0
	"PeakVal" em_doping
	"ValueAtDepth" base_doping
	"Depth" 1
	"Gauss" "Factor" 0.8)
(sdedr:define-refinement-window "RefEvalWin.AlBSF"
	"Rectangle"
	(position 0 0 (- THICKNESS))
	(position cellWidth cellLength (- THICKNESS) )
) ; window for the analytical profile placement.
(sdedr:define-analytical-profile-placement "PlaceAP.AlBSF"
	"Gauss.AlBSF"
	"RefEvalWin.AlBSF"
	"Positive"
	"NoReplace" "Eval"
) ; Place the analytical profile

; *** CONTACTS ***
(sdegeo:define-contact-set "em_contact" 	4 (color:rgb 1 0 0) "##") 	; em_contact
(sdegeo:define-contact-set "base_contact" 4 (color:rgb 0 0 1) "##") 	; base_contact

; Create the metal region
(sdegeo:create-cuboid (position 0 0 0) (position contactWidth cellLength 0.3) "Silver" "EmitterContact")
(sdegeo:delete-region (find-body-id (position (/ contactWidth 2) (/ cellLength 2) 0.15)))

; Activate contacts
(sdegeo:set-current-contact-set "em_contact")
(sdegeo:set-contact-faces (find-face-id (position (* contactWidth 0.5) (* cellLength 0.5) 0)) "em_contact")

(sdegeo:create-cuboid (position 0 0 (- THICKNESS) ) (position cellWidth cellLength (- (+ THICKNESS tBSF)) ) "Aluminum" "Region.BSF")
(sdegeo:delete-region (find-body-id (position (/ cellWidth 2) (/ cellLength 2) (- THICKNESS (- (/ tBSF 2))) )))
(sdegeo:set-current-contact-set "base_contact")
(sdegeo:set-contact-faces (find-face-id (position (* cellWidth 0.5) (* cellLength 0.5) (- THICKNESS))) "base_contact")

; *** OPTICAL GENERATION PROFILE ***
; NOTE: Make sure the window is defined in a direction normal to the optical generation profile!
(sdedr:define-refinement-window "opt_win"
	"Rectangle"
	(position contactWidth 0 0)
	(position cellWidth cellLength 0)
)
(sdedr:define-1d-external-profile "1d_opt_def" OptGenFile "Scale" 1.0 "Range" 0 THICKNESS "Erf" "Factor" 0) ;
(sdedr:define-analytical-profile-placement "opt_place" "1d_opt_def" "opt_win" "Negative" "NoReplace" "Eval")

; *** MESH ***
; * WHOLE DOMAIN
(sdedr:define-refeval-window "domain-ref" "Cuboid" (position 0 0 (- thicknessARC)) (position cellWidth cellLength (- THICKNESS)))
(sdedr:define-refinement-size "domain-ref-size" xmax ymax zmax xmin ymin zmin)
(sdedr:define-refinement-placement "domain-ref-pl" "domain-ref-size" "domain-ref")

; * p-n JUNCTION REFINEMENT * Used for an abrupt junction but shouldn't be necessary for a erf doping profile
(if (string=? dopingmodel "uniform")
	(begin
		(sdedr:define-refeval-window "junction-ref" "Cuboid" (position 0 0 (- (- emitterThickness 0.050)) ) (position cellWidth cellLength (-(+ emitterThickness 0.050))))
		(sdedr:define-refinement-size "junction-ref-size" (/ xmax 10) (/ ymax 10) (/ zmax 10) (/ xmin 10) (/ ymin 10) (/ zmin 10))
		(sdedr:define-refinement-placement "junction-ref-pl" "junction-ref-size" "junction-ref")
	)
)

; Doping refinement
(sdedr:define-refeval-window "RefWin.All"
	"Cuboid"
	(position 0 0 0)
	(position cellWidth cellLength (- THICKNESS) )
)
(sdedr:define-refinement-size "RefDef.All"
	(/ xmax 20) (/ ymax 20) (/ zmax 100)
	(/ xmin 20) (/ ymin 20) (/ zmin 100)
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

; Mesh refinement in the doped region (go up to 1 um)
(sdedr:define-refeval-window "junction-ref" "Cuboid" (position 0 0 0) (position cellWidth cellLength (- 1.2)))
(sdedr:define-refinement-size "junction-ref-size" (/ xmax 25) (/ ymax 25) (/ zmax 50) (/ xmin 50) (/ ymin 50) (/ zmin 100))
(sdedr:define-refinement-placement "junction-ref-pl" "junction-ref-size" "junction-ref")

;*********************** NOTE: add cond statement to define mesh refinement depending on the type of emitter
; Mesh refinement at the shunt
;(sdedr:define-refeval-window "shunt1_refine_window" "Rectangle" (position (- shunt1_pos_x1 0.05) 0 0) (position (+ shunt1_pos_x2 0.05) (+ ${dshunt_name} 0.05) 0) ) ; mesh refinement window slightly larger than the shunt region
;(if (> shuntDY 0)
;	(begin
;		(sdedr:define-refeval-window "shunt1_refine_window" "Cuboid" (position (- shunt1_pos_x1 0.05) (- shunt1_pos_y1 0.05) 0) (position (+ shunt1_pos_x2 0.05) (+ shunt1_pos_y2 0.05) (- (+ ${dshunt_name} 0.05))) ) ; mesh refinement window slightly larger than the shunt region
;		(sdedr:define-refinement-size "shunt-ref-size"
;			(/ shuntWidth 2) (/ shuntLength 2) (/ ${dshunt_name} 2)
;			(/ shuntWidth 10) (/ shuntLength 10) (/ ${dshunt_name} 10))
;	)
;)
;(sdedr:define-refinement-placement "shunt-ref-pl" "shunt-ref-size" "shunt1_refine_window")

(sdedr:define-refeval-window "RefWin.AlBSF"
	"Cuboid"
	(position 0 0 (- (- THICKNESS 1)))
	(position cellWidth cellLength (- THICKNESS))
)
(sdedr:define-refinement-size "RefDef.AlBSF"
	(/ xmax 10) (/ ymax 10) (/ zmax 20)
	(/ xmin 10) (/ ymin 10) (/ zmin 20)
)
(sdedr:define-refinement-placement "PlaceRF.AlBSF" "RefDef.AlBSF" "RefWin.AlBSF")

;(sdedr:define-refeval-window "RefWin.Optical"
;	"Cuboid"
;	(position 0 0 thicknessARC)
;	(position cellWidth cellLength (- (/ THICKNESS 4)))
;)
;(sdedr:define-refinement-size "RefDef.Optical"
;	(/ xmax 2) (/ ymax 2) (/ zmax 2)
;	(/ xmin 10) (/ ymin 10) (/ zmin 10)
;)
;(sdedr:define-refinement-placement "PlaceRF.Optical" "RefDef.Optical" "RefWin.Optical")

; SiNx refinement and SiNx/silicon refinement for interference calculation
;(sdedr:define-refeval-window "SiNx-Si-ref" "Rectangle" (position 0 (- thicknessARC) 0) (position cellWidth 0.050 0))
;(sdedr:define-refinement-size "SiNx-Si-ref-size" (/ xmax 500) (/ ymax 500) (/ xmin 10) (/ ymin 10))
;(sdedr:define-refinement-placement "SiNx-Si-ref-ref-pl" "SiNx-Si-ref-size" "SiNx-Si-ref-ref")

(display "  front contact refinement") (newline)
(sdedr:define-refinement-window "frontContact" "Cuboid"
	(position  0 0 0)
	(position  (+ contactWidth 2) cellLength (- THICKNESS))
)
(sdedr:define-refinement-size "frontContact"
	(/ contactWidth 2) (/ cellLength 10) (/ zmax 2)
	(/ contactWidth 10) (/ cellLength 20) (/ zmin 10)
)
(sdedr:define-refinement-placement "frontContact" "frontContact" "frontContact" )


; * BUILDING MESH
;(sde:build-mesh "snmesh" "-a -c boxmethod" "${nodenum}") ; save mesh files with name "nodnum"
(sde:build-mesh "snmesh" "-m 1000000 " "${nodenum}") ; save mesh files with name "nodnum"
