; *** Al-BSF solar cell for Na migration***
; Guillaume von Gastrow
; These command are written in Scheme programming language (based on LISP)

; *** INITIALIZATION ***
; clear structure
(sde:clear)

; New-replace-old option (default)
(sdegeo:set-default-boolean "ABA") ; subtract overlapping regions from the existing regions

; *** CHOICE OF DOPING PROFILE ***
(define dopingmodel "uniform")
;(define dopingmodel "erf")

; *** DEFINITIONS ***

; Initialize error file errfile.txt to no error (0)
(call-with-output-file "errfile.txt"
	(lambda (file)
		(write 0 file)))

; MATERIALS
(define CellMaterial "Silicon")
(define ARCMaterial "Si3N4")
(define shuntmat "Metal") ; Conductivity is added later as an external profile

; define contact length
(define Wcontact ${contact_length})
; define cell length
(define L ${device_length})
; wafer thickness
(define WSi 300)
; uniform emitter thickness (if selected)
(define Wem 0.6)
; The thickness of the Al BSF
(define tBSF 0.3)

(define err_flag 0)

; emitter region thickness (for constant emitter doping)
;(define Wem 100) ; Emitter depth 800 nm

; Silicon nitride thickness
(define tSiNx 75e-3)
; Shunt width
(define shuntw 0.01) ; um
; Shunt depths
(define ${dshunt_name} ${shunt_depth}) ; modified by the batch Sentaurus script based on the Na profile

; Position of shunt1
(define shunt1_pos_x1 (+ Wcontact (/ (- L Wcontact) 2)) ) ; left limit on the x-axis
(define shunt1_pos_x2 (+ (+ Wcontact (/ (- L Wcontact) 2) ) shuntw) ) ; right limit on the x-axis

(display shunt1_pos_x1)

; DOPING PARAMETERS
(define em_doping 1e19) ; Uniform doping value or surface dopant concentration, depending on doping model
(define base_doping 1e16) ; Silicon wafer base doping

; MESH PARAMETERS
(define xmax 10)
(define xmin 0.1)
(define ymax 1)
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
(sdegeo:create-rectangle (position Wcontact 0 0) (position L (- tSiNx) 0) ARCMaterial "SiNx-region")

; Choose between uniform doping and error function doping (ifelse statement, similar to switch syntax)
(cond
  ((string=? dopingmodel "uniform") ; condition 1, uniform doping
	(begin
		(newline)
		(display "Using an uniform doping profile")
		(newline)
		; create n-emitter
		(sdegeo:create-rectangle (position 0 0 0) (position L Wem 0) CellMaterial "Si-profile-region")
		; create p-emitter base
		(sdegeo:create-rectangle (position 0 Wem 0) (position L WSi 0) CellMaterial "Si-base-region")

		; *** DOPING ****
		; Emitter
		(sdedr:define-constant-profile "emitter-profile" "BoronActiveConcentration" em_doping)
		(sdedr:define-constant-profile-region "em-doping-placement" "emitter-profile" "Si-profile-region")

		; Base
		(sdedr:define-constant-profile "base-doping-profile" "PhosphorusActiveConcentration" base_doping)
		(sdedr:define-constant-profile-region "base-doping-placement" "base-doping-profile" "Si-base-region") ; place the base doping profile in the base region
		;(sdedr:define-constant-profile-placement "base-doping-placement" "base-doping-profile" "base-doping-window") ; not used because base-doping-window not defined
	)
  )
  ((string=? dopingmodel "erf") ; condition 2, error function doping
	(begin
		(display "Using an error function doping profile")
		(newline)
		; create Si region
		(sdegeo:create-rectangle (position 0 0 0) (position L WSi 0) CellMaterial "Si-profile-region") ; unique region on which an analytical doping profile will be added

		; *** DOPING ****
		; Emitter
		;(sdedr:define-constant-profile "emitter-profile" "PhosphorusActiveConcentration" em_doping) ; constant emitter doping
		; (sdedr:define-erf-profile "emitter-profile" "BoronActiveConcentration" "SymPos" 0.0 "MaxVal" (* 2 em_doping) "Length" 0.2 "erf" "factor" 1) ; error function doping profile. C0 is "MaxVal" (surface conc is C0/2) and "Length" is 2sqrt(Dt). See sde manual p. 222 and smesh manual p.114. Parameters based on profiles from Kerr et al in Fig. 4, JAP, 89, 2001.
		(sdedr:define-erf-profile "emitter-profile" "BoronActiveConcentration" "SymPos" 0.0 "MaxVal" (* 1 em_doping) "ValueAtDepth" base_doping Wem "erf" "factor" 1) ; error function doping profile. C0 is "MaxVal" (surface conc is C0/2) and "Length" is 2sqrt(Dt). See sde manual p. 222 and smesh manual p.114. Parameters based on profiles from Kerr et al in Fig. 4, JAP, 89, 2001.
		(sdedr:define-refinement-window "em-doping-window" "line" (position 0 0 0) (position L 0 0)) ; window for the analytical profile placement. In 2D, must be a line normal to the profile.
		(sdedr:define-analytical-profile-placement "em-doping-placement" "emitter-profile" "em-doping-window") ; Place the analytical profile

		; Base
		(sdedr:define-constant-profile "base-doping-profile" "PhosphorusActiveConcentration" base_doping)
		(sdedr:define-constant-profile-region "base-doping-placement" "base-doping-profile" "Si-profile-region") ; place the base doping profile in the base region
		;(sdedr:define-constant-profile-placement "base-doping-placement" "base-doping-profile" "base-doping-window") ; not used because base-doping-window not defined
	)
  )
)

; *** SHUNTS ***
; create shunts, ie conductive regions through the PN junction with recombination centers at the interface with silicon
; Later will need to place shunts with a do-loop across the width of the PN junction
; position of the shunt Wcontact+(L-Wcontact)/2

(display ${dshunt_name})
(display "\n")

; Define the shunts
(if (> shuntDY 0)
	(begin
		(for-each
			(lambda (SHUNTY INDEX)
				(begin
					(define REGION (string-append "shunt.Region." (number->string INDEX) ))
					(sdegeo:create-rectangle
						(position shunt1_pos_x1 SHUNTY 0)
						(position shunt1_pos_x2 (+ SHUNTY shuntDY) 0)
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

(sdegeo:create-rectangle (position 0 WSi 0) (position L (+ WSi tBSF) 0) "Aluminum" "Region.BSF")

; *** CONTACTS ***
; a) SET VERTICES
; 1st vertex on em_contact
(sdegeo:insert-vertex (position 0 0 0))
; 2nd vertex on em_contact (only part of the front side)
(sdegeo:insert-vertex (position Wcontact 0 0))
;(sdegeo:insert-vertex (position L 0 0))

; em_contact
(sdegeo:define-contact-set "em_contact" 4 (color:rgb 1 0 0) "##")
(sdegeo:set-current-contact-set "em_contact")
(sdegeo:define-2d-contact (find-edge-id (position (* Wcontact 0.5) 0 0)) "em_contact")
;(sdegeo:define-2d-contact (find-edge-id (position (* L 0.5) 0 0)) "em_contact")

; 1st vertex on base_contact
(sdegeo:insert-vertex (position 0 (+ WSi tBSF) 0))
; 2nd vertex on base_contact (contact covers the whole back side)
(sdegeo:insert-vertex (position L (+ WSi tBSF) 0))

; b) SET EDGE (DECLARATION, ACTIVATION AND DEFINITION)

; base_contact
(sdegeo:define-contact-set "base_contact" 4 (color:rgb 1 0 0) "##")
(sdegeo:set-current-contact-set "base_contact")
(sdegeo:define-2d-contact (find-edge-id (position (* L 0.5) (+ WSi tBSF) 0)) "base_contact")

; *** OPTICAL GENERATION PROFILE ***
; NOTE: Make sure the window is defined in a direction normal to the optical generation profile!
;(sdedr:define-refinement-window "opt_win" "Rectangle" (position (* Wcontact 0.8) 0 0) (position L (+ Wem Wbase) 0))
; (sdedr:define-refinement-window "opt_win" "Line" (position (* Wcontact 0.8) 0 0) (position L 0 0))
(sdedr:define-refinement-window "opt_win" "Line" (position (* Wcontact 1.0) 0 0) (position L 0 0))
;(sdedr:define-1d-external-profile "1d_opt_def" OptGenFile "Scale" 1.0 "Range" 0 1000 "Erf" "Factor" 0) ; Define the optical profile decreasing according to an error function
;(sdedr:define-1d-external-profile "1d_opt_def" OptGenFile "Scale" 1.0 "Range" 0 180 "Erf" "Length" 0) ; Define the optical profile decreasing according to an error function
; (sdedr:define-1d-external-profile "1d_opt_def" OptGenFile "Scale" 1.0 "Range" 0 WSi "Erf" "Length" 0) ; probably need to replace "Length" by "Factor" to avoid lateral spreading of the profile
(sdedr:define-1d-external-profile "1d_opt_def" OptGenFile "Scale" 1.0 "Range" 0 WSi "Erf" "Factor" 0) ;
; (sdedr:define-analytical-profile-placement "opt_place" "1d_opt_def" "opt_win" "Positive" "Replace" "Eval")
(sdedr:define-analytical-profile-placement "opt_place" "1d_opt_def" "opt_win" "Positive" "NoReplace" "Eval")

; *** MESH ***
; * WHOLE DOMAIN

(sdedr:define-refeval-window "domain-ref" "Rectangle" (position 0 (- tSiNx) 0) (position L WSi 0))
;(sdedr:define-refeval-window "domain-ref" "Rectangle" (position 0 0 0) (position L (+ Wem Wbase) 0))
(sdedr:define-refinement-size "domain-ref-size" xmax ymax xmin ymin)
(sdedr:define-refinement-placement "domain-ref-pl" "domain-ref-size" "domain-ref")



; * p-n JUNCTION REFINEMENT * Used for an abrupt junction but shouldn't be necessary for a erf doping profile
(if (string=? dopingmodel "uniform")
	(begin
		(sdedr:define-refeval-window "junction-ref" "Rectangle" (position 0 (- Wem 0.050) 0) (position L (+ Wem 0.050) 0))
		(sdedr:define-refinement-size "junction-ref-size" (/ xmax 10) (/ ymax 10) (/ xmin 10) (/ ymin 10))
		(sdedr:define-refinement-placement "junction-ref-pl" "junction-ref-size" "junction-ref")
	)
)

; Mesh refinement in the doped region (go up to 1 um)
(sdedr:define-refeval-window "junction-ref" "Rectangle" (position 0 0 0) (position L 1.2 0))
(sdedr:define-refinement-size "junction-ref-size" (/ xmax 100) (/ ymax 100) (/ xmin 100) (/ ymin 100))
(sdedr:define-refinement-placement "junction-ref-pl" "junction-ref-size" "junction-ref")

;*********************** NOTE: add cond statement to define mesh refinement depending on the type of emitter
; Mesh refinement at the shunt
;(sdedr:define-refeval-window "shunt1_refine_window" "Rectangle" (position (- shunt1_pos_x1 0.05) 0 0) (position (+ shunt1_pos_x2 0.05) (+ ${dshunt_name} 0.05) 0) ) ; mesh refinement window slightly larger than the shunt region
(if (> shuntDY 0)
	(begin
		(sdedr:define-refeval-window "shunt1_refine_window" "Rectangle" (position (- shunt1_pos_x1 0.05) (- tSiNx) 0) (position (+ shunt1_pos_x2 0.05) (+ ${dshunt_name} 0.05) 0) ) ; mesh refinement window slightly larger than the shunt region
		;(sdedr:define-refinement-size "shunt-ref-size" (/ xmax 20) (/ ymax 20) (/ xmin 20) (/ ymin 20))
		(sdedr:define-refinement-size "shunt-ref-size" (/ shuntw 2) (/ ${dshunt_name} 2) (/ shuntw 20) (/ ${dshunt_name} 20))
	)
)


(sdedr:define-refinement-placement "shunt-ref-pl" "shunt-ref-size" "shunt1_refine_window")

; SiNx refinement and SiNx/silicon refinement for interference calculation
;(sdedr:define-refeval-window "SiNx-Si-ref" "Rectangle" (position 0 (- tSiNx) 0) (position L 0.050 0))
;(sdedr:define-refinement-size "SiNx-Si-ref-size" (/ xmax 500) (/ ymax 500) (/ xmin 10) (/ ymin 10))
;(sdedr:define-refinement-placement "SiNx-Si-ref-ref-pl" "SiNx-Si-ref-size" "SiNx-Si-ref-ref")

(display "  front contact refinement") (newline)
(sdedr:define-refinement-window "frontContact" "Rectangle"
	(position  0 0 0)
	(position  (+ Wcontact 2) 1 0)
)
(sdedr:define-refinement-size "frontContact"
	(/ Wcontact 10) (/ Wcontact 50)
	(/ Wcontact 10) (/ Wcontact 50)
)
(sdedr:define-refinement-placement "frontContact" "frontContact" "frontContact" )


; * BUILDING MESH
(sde:build-mesh "snmesh" "-a -c boxmethod" "${nodenum}") ; save mesh files with name "nodnum"
