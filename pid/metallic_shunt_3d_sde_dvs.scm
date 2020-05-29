; *** Al-BSF solar cell for Na migration***
; Erick Martinez Loran
; These command are written in Scheme programming language (based on LISP)

; *** INITIALIZATION ***
; clear structure
(sde:clear)

; MATERIALS
(define CellMaterial 			"Silicon")
(define ARCMaterial 			"Si3N4")
(define shuntmat 					"Metal")
(define contactWidth 			${contact_length})
(define cellWidth 				${cell_width})
(define cellLength 				${cell_length})
(define THICKNESS 				200)
(define emitterThickness 	0.6)
(define thicknessARC 			75e-3)

; The thickness of the Al BSF
(define tBSF 			0.3)

(define SHUNT_DEPTH 			${shunt_depth}) ; modified by the batch Sentaurus script based on the Na profile


; New-replace-old option (default)
(sdegeo:set-default-boolean "ABA") ; subtract overlapping regions from the existing regions

; DOPING PARAMETERS
(define em_doping 1e19) ; Uniform doping value or surface dopant concentration, depending on doping model
(define base_doping 1e16) ; Silicon wafer base doping

; MESH PARAMETERS
(define xmax 10)
(define xmin 0.1)
(define ymax 10)
(define ymin 0.1)
(define zmax 20)
(define zmin 0.1)
; Import optical generation file
(define OptGenFile "${optical_generation_file}" )

; *** SHUNTS ***
; create shunts, ie conductive regions through the PN junction with recombination centers at the interface with silicon
${shunts}

; create SiNx film (leave space for the metal contact before)
(sdegeo:create-cuboid (position (/ contactWidth 2) 0 0) (position (- cellWidth (/ contactWidth 2)) cellLength thicknessARC) ARCMaterial "SiNx-region")

(if (> SHUNT_DEPTH 0)
	(begin
		; Old-replace-new option (default)
		(sdegeo:set-default-boolean "BAB") ; subtract overlapping regions from the existing regions
	)
)


; *** GEOMETRY ***
; convention: x=length y=thickness
; *** Silicon ***
; create Si region
(sdegeo:create-cuboid
	(position 0 0 0)
	(position cellWidth cellLength (- THICKNESS) )
	CellMaterial "Si-profile-region"
)


(if (> SHUNT_DEPTH 0)
	(begin
		; New-replace-old option (default)
		(sdegeo:set-default-boolean "ABA") ; subtract overlapping regions from the existing regions
	)
)




; *** DOPING ****
; Base
(sdedr:define-constant-profile
	"base-doping-profile"
	"BoronActiveConcentration"
	base_doping
)
(sdedr:define-constant-profile-region
	"base-doping-placement"
	"base-doping-profile"
	"Si-profile-region"
) ; place the base doping profile in the base region
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

; Aluminum Back Surface Field
(sdedr:define-gaussian-profile "Gauss.AlBSF"
	"AluminumActiveConcentration"
	"PeakPos" 0
	"PeakVal" em_doping
	"ValueAtDepth" base_doping
	"Depth" 1
	"Gauss" "Factor" 0.5)
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
(sdegeo:create-cuboid (position 0 0 0) (position (/ contactWidth 2) cellLength 1) "Silver" "EmitterContact")
(sdegeo:create-cuboid (position (- cellWidth (/ contactWidth 2)) 0 0) (position cellWidth cellLength 1) "Silver" "EmitterContact")
(sdegeo:delete-region (find-body-id (position (/ contactWidth 4) (/ cellLength 2) 0.5)))
(sdegeo:delete-region (find-body-id (position (- cellWidth (/ contactWidth 4)) (/ cellLength 2) 0.5)))

; Activate contacts
(sdegeo:set-current-contact-set "em_contact")
(sdegeo:set-contact-faces (find-face-id (position (* contactWidth 0.25) (* cellLength 0.5) 0)) "em_contact")
(sdegeo:set-contact-faces (find-face-id (position (- cellWidth (* contactWidth 0.25)) (* cellLength 0.5) 0)) "em_contact")

(sdegeo:create-cuboid (position 0 0 (- THICKNESS) ) (position cellWidth cellLength (- (+ THICKNESS tBSF)) ) "Aluminum" "Region.BSF")
(sdegeo:delete-region (find-body-id (position (/ cellWidth 2) (/ cellLength 2) (- (+ THICKNESS (/ tBSF 2))) )))
(sdegeo:set-current-contact-set "base_contact")
(sdegeo:set-contact-faces (find-face-id (position (* cellWidth 0.5) (* cellLength 0.5) (- THICKNESS))) "base_contact")

; *** MESH ***
; * WHOLE DOMAIN
; Doping refinement
(sdedr:define-refeval-window "RefWin.All"
	"Cuboid"
	(position 0 0 0)
	(position cellWidth cellLength (- THICKNESS) )
)
(sdedr:define-refinement-size "RefDef.All"
	xmax ymax zmax
	xmin ymin zmin
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

; Mesh refinement around the junction (go up down to twice the shunt depth)
(sdedr:define-refeval-window "RefWin.Junction" "Cuboid" (position 0 0 0) (position cellWidth cellLength 5))
(sdedr:define-refinement-size "RefDef.Junction" (/ xmax 10) (/ ymax 10) (/ ymax 10) (/ xmin 10) (/ ymin 10) (/ zmin 10))
(sdedr:define-refinement-placement "PlacementRF.Junction" "RefDef.Junction" "RefWin.Junction")

;*********************** NOTE: add cond statement to define mesh refinement depending on the type of emitter
; Mesh refinement around the shunt
${shunt_refinement}


(sdedr:define-refeval-window "RefWin.AlBSF"
	"Cuboid"
	(position 0 0 (- (- THICKNESS 1)))
	(position cellWidth cellLength (- THICKNESS))
)
(sdedr:define-refinement-size "RefDef.AlBSF"
	(/ xmax 10) (/ ymax 10) (/ zmax 10)
	(/ xmin 10) (/ ymin 10) (/ zmin 10)
)
(sdedr:define-refinement-placement "PlaceRF.AlBSF" "RefDef.AlBSF" "RefWin.AlBSF")

; *** OPTICAL GENERATION PROFILE ***
; NOTE: Make sure the window is defined in a direction normal to the optical generation profile!
(sdedr:define-refinement-window "RefWin.Optical"
	"Rectangle"
	(position (/ contactWidth 2) 0 0)
	(position (- cellWidth (/ contactWidth 2)) cellLength 0)
)
(sdedr:define-1d-external-profile "RefDef.Optical_1D" OptGenFile "Scale" 1.0 "Range" 0 THICKNESS "Erf" "Factor" 0) ;
(sdedr:define-analytical-profile-placement "PlaceRF.Optical_1D" "RefDef.Optical_1D" "RefWin.Optical" "Negative" "NoReplace" "Eval")


(display "  front contact refinement") (newline)
(sdedr:define-refinement-window "RefWin.frontContact1" "Cuboid"
	(position  0 0 0)
	(position  (+ contactWidth 1) cellLength (- 10))
)

(sdedr:define-refinement-window "RefWin.frontContact2" "Cuboid"
	(position  (- cellWidth (+ contactWidth 1)) 0 0)
	(position  cellWidth cellLength (- 10))
)

(sdedr:define-refinement-size "RefDef.frontContact"
	(/ contactWidth 10) (/ cellLength 10) (/ zmax 10)
	(/ contactWidth 50) (/ cellLength 50) (/ zmin 10)
)
(sdedr:define-refinement-placement "PlaceRef.frontContact1" "RefDef.frontContact" "RefWin.frontContact1" )
(sdedr:define-refinement-placement "PlaceRef.frontContact2" "RefDef.frontContact" "RefWin.frontContact2" )

(display "Finished refinements. Ready to mesh!") (newline)

; * BUILDING MESH
;(sde:build-mesh "snmesh" "-a -c boxmethod" "${nodenum}") ; save mesh files with name "nodnum"
(sde:build-mesh " -offset -m 1000000 --threads 40 " "${nodenum}") ; save mesh files with name "nodnum"
