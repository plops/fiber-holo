ccl::*shared-libraries*
ccl::*eeps* ;; hash table with external functions

(progn (use-interface-dir :andor)
	(open-shared-library "libandor-emu.so"))
(progn
 
 (defun get-available-cameras ()
   (rlet ((n :int))
     (assert (= #$DRV_SUCCESS
		(#_GetAvailableCameras
		 (%int-to-ptr (%address-of n)))))
     (%get-signed-long (%int-to-ptr (%address-of n)))))
 (defun get-camera-handle (&optional (i 0))
   (rlet ((n :int))
     (assert (= #$DRV_SUCCESS
		(#_GetCameraHandle i (%int-to-ptr (%address-of n)))))
     (%get-signed-long (%int-to-ptr (%address-of n)))))
 (defun set-current-camera (handle)
   (#_SetCurrentCamera handle))
 (defun initialize (&optional (dir "/usr/local/etc/andor/"))
   (with-cstrs ((cdir dir))
     (assert (= #$DRV_SUCCESS (#_Initialize cdir)))))
 (defun start-acquisition ()
   (assert (= #$DRV_SUCCESS (#_StartAcquisition))))
 (defun get-status ()
   (rlet ((status :int))
     (assert (= #$DRV_SUCCESS
		(#_GetStatus
		 (%int-to-ptr (%address-of status)))))
     (%get-signed-long (%int-to-ptr (%address-of status)))))
 (defun get-acquired-data16 ()
   (let* ((w 512)
	  (h 512)
	  (n (* w h)))
     (multiple-value-bind (a ap) (make-heap-ivector n ;:unsigned-halfword 
						    '(unsigned-byte 16)
						    )
       (let ((a (#_GetAcquiredData16 ap n)))
	 (unless (= a #$DRV_SUCCESS)
	   (break "GetAcquiredData16 returned ~a." a)))
       (prog1 
	   (let* ((b1 (make-array n :element-type '(unsigned-byte 16)))
		  (b (make-array (list h w) :element-type '(unsigned-byte 16)
				 :displaced-to b1)))
	     (dotimes (i n)
	       (setf (aref b1 i) (aref a i)))
	     b)
	 (dispose-heap-ivector a))))))
#+nil
(progn
  (initialize)
  (start-acquisition)
  (external-call "print_cam" :void)
  (get-status)
  (external-call "print_cam" :void)
  (defparameter *bla*   (get-acquired-data16))
    (external-call "print_cam" :void))

#+nil
(external-call "print_cam" :void)

(load "/home/martin/quicklisp/setup.lisp")

(progn
  (ql:quickload :cl-glut)
  (ql:quickload :cl-opengl)
  (ql:quickload :cl-glu))

(defpackage :g (:use :cl :gl :ccl))

(in-package :g)

(let ((rot 0))
 (defun draw ()
   (clear :color-buffer :depth-buffer)
   (color 1 1 1)
   (incf rot 1)
   (if (< 360 rot)
       (setf rot 0))
   (with-pushed-matrix
     (translate 0 0 -5)
     (rotate rot 1 0 0)
     (line-width 4)
     (let ((r 1.5))
       (glut:wire-sphere r 15 7)
       (color 0 0 0 1)
       (glut:solid-sphere (* .99 r) 15 7)))

   (with-pushed-matrix
     (translate (sin (* 8 2 pi (/ rot 360))) 0 -2.2)
     (color 1 1 1)
     (rect .14 -3 -.14 3))
   (with-pushed-matrix
     (translate -.5 -.5 -1.4)
     (color 1 1 1)
     (let ((obj (first (gen-textures 1))))
       (bind-texture :texture-2d obj)
       
       (tex-parameter :texture-2d :texture-min-filter :nearest)
       (tex-parameter :texture-2d :texture-mag-filter :nearest)
       (start-acquisition) 
       (let ((b1 (make-array (* 512 512) :element-type '(unsigned-byte 16)
			     :displaced-to (get-acquired-data16))))
	 (multiple-value-bind (a ap)
	     (make-heap-ivector (* 512 512)
				'(unsigned-byte 16))
	   (dotimes (i (* 512 512))
	     (setf (aref a i) (* 500 (aref b1 i))))
	   (tex-image-2d :texture-2d 0 :rgba 512 512 0 :green :unsigned-short
			 ap) 
	   (dispose-heap-ivector a)))
       (enable :texture-2d)
       (with-primitive :quads
	 (vertex 0 0) (tex-coord 0 0)
	 (vertex 0 1) (tex-coord 0 1)
	 (vertex 1 1) (tex-coord 1 1)
	 (vertex 1 0) (tex-coord 1 0))
       (disable :texture-2d)
       (delete-textures (list obj))))
   (flush)
 ;;  (sleep (/ 61.3))
;   (cl-opengl-bindings:wait-sync )
   (glut:post-redisplay)))


(progn 

  (defclass clip-window (glut:window)
    ()
    (:default-initargs :pos-x 853 :pos-y 0 :width 512 :height 768
		       :mode '(:single :rgb :depth) :title "clip.lisp"))

  (defmethod glut:display-window :before ((w clip-window))
    (clear-color 0 0 0 0)
    (enable :depth-test)
    (shade-model :flat))

  (defmethod glut:display ((w clip-window))
    (draw))

  (defmethod glut:reshape ((w clip-window) width height)
    (viewport 0 0 width height)
    (matrix-mode :projection)
    (load-identity)
    (glu:perspective 60 (/ width height) 1 20)
    (matrix-mode :modelview))

  (defmethod glut:keyboard ((w clip-window) key x y)
    (declare (ignore x y))
    (when (eql key #\Esc)
      (glut:destroy-current-window)))

  (defun rb-cam ()
    (glut:display-window (make-instance 'clip-window)))
  )

#+nil
(ccl:process-run-function "gl-display" #'rb-cam)