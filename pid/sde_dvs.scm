; *** Al-BSF solar cell for Na migration***
; Erick Martinez Loran & Guillaume von Gastrow

; *** INITIALIZATION ***
; clear structure
(sde:clear)

; New-replace-old option (default)
; subtract overlapping regions from the existing regions
(sdegeo:set-default-boolean "ABA") 

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
(define shuntmat "Silicon") ; Conductivity is added later as an external profile

; define contact length
(define Wcontact ${contact_length})
; define cell length
(define Wdevice ${device_length})
; wafer thickness
(define tSi 300)
; uniform emitter thickness (if selected)
(define tEmitter 0.6)
; Error flag
(define err_flag 0)
; Silicon nitride thickness
(define tSiNx 75e-3)
; Shunt width
(define shuntw 0.01) ; um
; Shunt depths
; modified by the batch Sentaurus script based on the Na profile
(define ${dshunt_name} ${shunt_depth})
; Position of shunt1
; left limit on the x-axis
(define shunt1_pos_x1 (+ Wcontact (/ (- Wdevice Wcontact) 2)) ) 
; right limit on the x-axis
(define shunt1_pos_x2 (+ (+ Wcontact (/ (- Wdevice Wcontact) 2) ) shuntw) ) 

; DOPING PARAMETERS
; Uniform doping value or surface dopant concentration, depending on doping model
(define em_doping 1e19) 
; Silicon wafer base doping
(define base_doping 1e16) 

; MESH PARAMETERS
(define xmax 10)
(define xmin 0.1)
(define ymax 1)
(define ymin 0.1)
; Import optical generation file
(define OptGenFile "${optical_generation_file}" )
; Import shunt conductivity file
(define shunt_profile "${shunt_profile_file}")
(display shunt_profile)
(display shunt1_pos_x1)

; *** GEOMETRY ***
; convention: x=length y=thickness

; create SiNx film (leave space for the metal contact before)
(sdegeo:create-rectangle (position Wcontact 0 0) (position Wdevice (- tSiNx) 0) ARCMaterial "SiNx-region")

; Choose between uniform doping and error function doping (ifelse statement, similar to switch syntax)
(cond
  ((string=? dopingmodel "uniform") ; condition 1, uniform doping
	(begin
		(newline)
		(display "Using an uniform doping profile")
		(newline)
		; create n-emitter
		(sdegeo:create-rectangle (position 0 0 0) (position Wdevice tEmitter 0) CellMaterial "Region.Emitter")
		; create p-emitter base
		(sdegeo:create-rectangle (position 0 tEmitter 0) (position Wdevice tSi 0) CellMaterial "Region.Base")

		; *** DOPING ****
		; Emitter
		(sdedr:define-constant-profile "Prof.Emitter" "BoronActiveConcentration" em_doping)
		(sdedr:define-constant-profile-region "Placement.Emitter" "Prof.Emitter" "Region.Emitter")

		; Base
		(sdedr:define-constant-profile "Prof.Base" "PhosphorusActiveConcentration" base_doping)
		(sdedr:define-constant-profile-region "Placement.Base" "Prof.Base" "Region.Base") ; place the base doping profile in the base region
	)
  )
  ((string=? dopingmodel "erf") ; condition 2, error function doping
	(begin
		(display "Using an error function doping profile")
		(newline)
		; create Si region
		(sdegeo:create-rectangle (position 0 0 0) (position Wdevice tSi 0) CellMaterial "Region.Emitter") ; unique region on which an analytical doping profile will be added

		; *** DOPING ****
		; Emitter
		; (sdedr:define-erf-profile "Prof.Emitter" "BoronActiveConcentration" "SymPos" 0.0 "MaxVal" (* 2 em_doping) "Length" 0.2 "erf" "factor" 1) ; error function doping profile. C0 is "MaxVal" (surface conc is C0/2) and "Length" is 2sqrt(Dt). See sde manual p. 222 and smesh manual p.114. Parameters based on profiles from Kerr et al in Fig. 4, JAP, 89, 2001.
		(sdedr:define-erf-profile "Prof.Emitter" "BoronActiveConcentration" "SymPos" 0.0 "MaxVal" (* 1 em_doping) "ValueAtDepth" base_doping tEmitter "erf" "factor" 1) ; error function doping profile. C0 is "MaxVal" (surface conc is C0/2) and "Length" is 2sqrt(Dt). See sde manual p. 222 and smesh manual p.114. Parameters based on profiles from Kerr et al in Fig. 4, JAP, 89, 2001.
		(sdedr:define-refinement-window "BaseLine.Emitter" "line" (position 0 0 0) (position Wdevice 0 0)) ; window for the analytical profile placement. In 2D, must be a line normal to the profile.
		(sdedr:define-analytical-profile-placement "Placement.Emitter" "Prof.Emitter" "BaseLine.Emitter") ; Place the analytical profile

		; Base
		(sdedr:define-constant-profile "Prof.Base" "PhosphorusActiveConcentration" base_doping)
		(sdedr:define-constant-profile-region "Placement.Base" "Prof.Base" "Region.Emitter") ; place the base doping profile in the base region
	)
  )
)

; *** SHUNTS ***
; create shunts, ie conductive regions through the PN junction with recombination centers at the interface with silicon
; Later will need to place shunts with a do-loop across the width of the PN junction
; position of the shunt Wcontact+(Wdevice-Wcontact)/2
(display ${dshunt_name})
(newline)

; Handle error if the shunt depth is 0 um by sending an error flag through an external file
(if (= ${dshunt_name} 0)
	(begin
		; Write value 1 into a file errfile.txt to indicate error
		(call-with-output-file "errfile.txt"
			(lambda (file)
				(write 1 file)
			)
		)
	)
)

; Define region for the shunt
;(sdegeo:create-rectangle (position (+ Wcontact(/(- Wdevice Wcontact)2)) 0 0) (position (+ (+ Wcontact(/(- Wdevice Wcontact)2)) shuntw) ${dshunt_name} 0) shuntmat "Region.Shunt1")
(sdegeo:create-rectangle
	(position shunt1_pos_x1 0 0)
	(position shunt1_pos_x2 ${dshunt_name} 0)
	shuntmat
	"Region.Shunt1") ; Check if this command is really needed, since the window definition overwrites this region (maybe only useful for plotting).
(sdedr:define-constant-profile "Prof.Shunt" "BoronActiveConcentration" em_doping)
(sdedr:define-constant-profile-region "Placement.ShuntDoping" "Prof.Emitter" "Region.Shunt1")
; Define external conductivity profile for the shunt
;(sdedr:define-1d-external-profile "shunt_cond" shunt_profile "Scale" 1.0 "Erf" "Length" 0)
(sdedr:define-1d-external-profile
	"DeepLevels"
	shunt_profile
	"Range" 0 ${dshunt_name}
	"Erf" "Length" 0.0)
;(sdedr:define-1d-external-profile "shunt_cond" shunt_profile "Scale" 1.0 "Erf" "Factor" 0)
; Define window for shunt ("Line" should be used, sde manual p. 561 + the window must be normal to the profile di0rection)
; (sdedr:define-refeval-window "Region.Shunt1" "Line" (position shunt1_pos_x1 0 0) (position shunt1_pos_x2 0 0) )
(sdedr:define-refeval-window
	"Region.Shunt1"
	"Line"
	(position shunt1_pos_x1 0 0)
	(position shunt1_pos_x2 0 0) )
; Define conductivity placement
(sdedr:define-analytical-profile-placement
	"shunt1_placement"
	"DeepLevels"
	"Region.Shunt1"
	"Positive"
	"Replace"
	"Eval" "" 0.0
	"evalwin")

;(sdegeo:create-rectangle (position (- Wdevice 0.1) 0 0) (position (+ (- Wdevice 0.1) shuntw) (* 0.9 tEmitter) 0) shuntmat "shunt-region")
;(sdegeo:create-rectangle (position (+ Wcontact(/(- Wdevice Wcontact)2) 10) 0 0) (position (+ (+ Wcontact(/(- Wdevice Wcontact)2) 10) shuntw) (* 1.1 tEmitter) 0) shuntmat "shunt-region2")
;(sdegeo:create-rectangle (position (+ Wcontact(/(- Wdevice Wcontact)2) 15) 0 0) (position (+ (+ Wcontact(/(- Wdevice Wcontact)2) 15) shuntw) (* 1.5 tEmitter) 0) shuntmat "shunt-region3")
;(sdegeo:create-rectangle (position (- (+ Wcontact(/(- Wdevice Wcontact)2)) 12) 0 0) (position (-(+ (+ Wcontact(/(- Wdevice Wcontact)2)) shuntw)12) (* 1.3 tEmitter) 0) shuntmat "shunt-region4")

; *** CONTACTS ***
; a) SET VERTICES
; 1st vertex on em_contact
(sdegeo:insert-vertex (position 0 0 0))
; 2nd vertex on em_contact (only part of the front side)
(sdegeo:insert-vertex (position Wcontact 0 0))
;(sdegeo:insert-vertex (position Wdevice 0 0))

; em_contact
(sdegeo:define-contact-set "em_contact" 4 (color:rgb 1 0 0) "##")
(sdegeo:set-current-contact-set "em_contact")
(sdegeo:define-2d-contact (find-edge-id (position (* Wcontact 0.5) 0 0)) "em_contact")
;(sdegeo:define-2d-contact (find-edge-id (position (* Wdevice 0.5) 0 0)) "em_contact")

; 1st vertex on base_contact
(sdegeo:insert-vertex (position 0 tSi 0))
; 2nd vertex on base_contact (contact covers the whole back side)
(sdegeo:insert-vertex (position Wdevice tSi 0))

; b) SET EDGE (DECLARATION, ACTIVATION AND DEFINITION)

; base_contact
(sdegeo:define-contact-set "base_contact" 4 (color:rgb 1 0 0) "##")
(sdegeo:set-current-contact-set "base_contact")
(sdegeo:define-2d-contact (find-edge-id (position (* Wdevice 0.5) tSi 0)) "base_contact")

; *** OPTICAL GENERATION PROFILE ***
; NOTE: Make sure the window is defined in a direction normal to the optical generation profile!
;(sdedr:define-refinement-window "opt_win" "Rectangle" (position (* Wcontact 0.8) 0 0) (position Wdevice (+ tEmitter Wbase) 0))
; (sdedr:define-refinement-window "opt_win" "Line" (position (* Wcontact 0.8) 0 0) (position Wdevice 0 0))
(sdedr:define-refinement-window "opt_win" "Line" (position (* Wcontact 1.0) 0 0) (position Wdevice 0 0))
;(sdedr:define-1d-external-profile "1d_opt_def" OptGenFile "Scale" 1.0 "Range" 0 1000 "Erf" "Factor" 0) ; Define the optical profile decreasing according to an error function
;(sdedr:define-1d-external-profile "1d_opt_def" OptGenFile "Scale" 1.0 "Range" 0 180 "Erf" "Length" 0) ; Define the optical profile decreasing according to an error function
; (sdedr:define-1d-external-profile "1d_opt_def" OptGenFile "Scale" 1.0 "Range" 0 tSi "Erf" "Length" 0) ; probably need to replace "Length" by "Factor" to avoid lateral spreading of the profile
(sdedr:define-1d-external-profile "1d_opt_def" OptGenFile "Scale" 1.0 "Range" 0 tSi "Erf" "Factor" 0) ;
; (sdedr:define-analytical-profile-placement "opt_place" "1d_opt_def" "opt_win" "Positive" "Replace" "Eval")
(sdedr:define-analytical-profile-placement "opt_place" "1d_opt_def" "opt_win" "Positive" "NoReplace" "Eval")

; *** MESH ***
; * WHOLE DOMAIN

(sdedr:define-refeval-window "domain-ref" "Rectangle" (position 0 (- tSiNx) 0) (position Wdevice tSi 0))
;(sdedr:define-refeval-window "domain-ref" "Rectangle" (position 0 0 0) (position Wdevice (+ tEmitter Wbase) 0))
(sdedr:define-refinement-size "domain-ref-size" xmax ymax xmin ymin)
(sdedr:define-refinement-placement "domain-ref-pl" "domain-ref-size" "domain-ref")



; * p-n JUNCTION REFINEMENT * Used for an abrupt junction but shouldn't be necessary for a erf doping profile
(if (string=? dopingmodel "uniform")
	(begin
		(sdedr:define-refeval-window "junction-ref" "Rectangle" (position 0 (- tEmitter 0.050) 0) (position Wdevice (+ tEmitter 0.050) 0))
		(sdedr:define-refinement-size "junction-ref-size" (/ xmax 10) (/ ymax 10) (/ xmin 10) (/ ymin 10))
		(sdedr:define-refinement-placement "junction-ref-pl" "junction-ref-size" "junction-ref")
	)
)

; Mesh refinement in the doped region (go up to 1 um)
(sdedr:define-refeval-window "junction-ref" "Rectangle" (position 0 0 0) (position Wdevice 1.2 0))
(sdedr:define-refinement-size "junction-ref-size" (/ xmax 100) (/ ymax 100) (/ xmin 100) (/ ymin 100))
(sdedr:define-refinement-placement "junction-ref-pl" "junction-ref-size" "junction-ref")

;*********************** NOTE: add cond statement to define mesh refinement depending on the type of emitter
; Mesh refinement at the shunt
;(sdedr:define-refeval-window "shunt1_refine_window" "Rectangle" (position (- shunt1_pos_x1 0.05) 0 0) (position (+ shunt1_pos_x2 0.05) (+ ${dshunt_name} 0.05) 0) ) ; mesh refinement window slightly larger than the shunt region
(sdedr:define-refeval-window "shunt1_refine_window" "Rectangle" (position (- shunt1_pos_x1 0.05) (- tSiNx) 0) (position (+ shunt1_pos_x2 0.05) (+ ${dshunt_name} 0.05) 0) ) ; mesh refinement window slightly larger than the shunt region
;(sdedr:define-refinement-size "shunt-ref-size" (/ xmax 20) (/ ymax 20) (/ xmin 20) (/ ymin 20))
(sdedr:define-refinement-size "shunt-ref-size" (/ shuntw 2) (/ ${dshunt_name} 2) (/ shuntw 20) (/ ${dshunt_name} 20))

(sdedr:define-refinement-placement "shunt-ref-pl" "shunt-ref-size" "shunt1_refine_window")

; SiNx refinement and SiNx/silicon refinement for interference calculation
;(sdedr:define-refeval-window "SiNx-Si-ref" "Rectangle" (position 0 (- tSiNx) 0) (position Wdevice 0.050 0))
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
