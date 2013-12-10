(load "/home/martin/quicklisp/setup.lisp")
(progn
  (ql:quickload :cl-glut)
  (ql:quickload :cl-opengl)
  (ql:quickload :cl-glu))

(defpackage :g (:use :cl :gl :ccl))
(in-package :g)

#.(load "/home/martin/src/ccl/library/serial-streams.lisp")

(defvar *serial* nil)
#+nil
(defparameter *serial*
  (ccl::make-serial-stream "/dev/ttyACM0"
                           ;:format 'character
                           :baud-rate 9600
                           :parity nil
                           :char-bits 8
                           :stop-bits 1
                           :flow-control nil))

#+nil
(close *serial*) 

;uvcdynctrl -d /dev/video1  -s "Exposure, Auto" 1
;uvcdynctrl -d /dev/video1 -s "Exposure (Absolute)" 900

(defun write-pgm (fn a)
  (destructuring-bind (h w) (array-dimensions a) 
    (with-open-file (s fn :direction :output
                       :if-exists  :supersede
                       :if-does-not-exist :create)
      (format s "P5~%~a ~a~%255~%" w h))
    (with-open-file (s fn :direction :output
                       :if-exists  :append
                       :if-does-not-exist :create
                       :element-type '(unsigned-byte 8))
      (dotimes (i w)
        (dotimes (j h)
         (write-byte (aref a j i) s))))))

#+nil
(progn
  (format *serial* "v-200~%" voltage)
  (force-output *serial*))


(defun write-mirror (voltage)
  (format *serial* "v~d~%" voltage)
  (force-output *serial*)
  (sleep .1)
  (list 
   (read-line *serial*)
   (read-line *serial*)))

#+nil
(write-mirror -35)

#.(progn
  (use-interface-dir :fftw3)
  (open-shared-library "libfftw3.so"))

#+nil
(let ((border .7))
 (loop for i below 20 collect 
      (let ((r (/ i 20d0)))
	)
      ))

(defparameter *window* 
  (let ((a (make-array (list 512 512) :element-type 'double-float))
	(border .5))
    (dotimes (i 512)
      (dotimes (j 512)
	(setf (aref a j i)
	      (*  (expt -1d0 (+ i j))
	       (let* ((x (/ (- i 256) 256d0))
		      (y (/ (- j 256) 256d0))
		      (r (sqrt (+ (* x x) (* y y)))
			))
		 (if (<  r (- 1 border)) 
		     1d0
		     (let ((rr (* (/ border) (- r (- 1 border)))))
		       (* .5 (+ 1 (cos (* pi rr)))))))))))
    a)) 

(defun ft (a)
   (destructuring-bind (h w) (array-dimensions a)
     (multiple-value-bind (ind indp) (ccl:make-heap-ivector 2 '(signed-byte 32))
       (setf (aref ind 0) h
	     (aref ind 1) w)
       (multiple-value-bind (b bp) (ccl:make-heap-ivector (* 2 h w) 'double-float)
	 (let ((a1 (make-array (* w h) :element-type (array-element-type a)
			       :displaced-to a))
	       (w1 (make-array (* 512 512) :element-type (array-element-type *window*)
			       :displaced-to *window*)))
	   (if (complexp (aref a1 0))
	       (dotimes (i (length a1))
		 (setf (aref b (* 2 i))       (* 1d0 (realpart (aref a1 i)))
		       (aref b (+ 1 (* 2 i))) (* 1d0 (imagpart (aref a1 i)))))
	       (dotimes (i (length a1))
		 (setf (aref b (* 2 i))       (* (aref w1 i) (aref a1 i))
		       (aref b (+ 1 (* 2 i))) 0d0))))
	 (multiple-value-bind (c cp) (ccl:make-heap-ivector (* 2 h w) 'double-float)	
	   (prog1
	       (progn
		 (let ((plan (#_fftw_plan_dft 2 indp bp cp #$FFTW_FORWARD #$FFTW_ESTIMATE)))
		   (when (%null-ptr-p plan)
		     (error "fftw_plan_dft ~a returned null." (list bp cp)))
		   (unwind-protect
			(#_fftw_execute plan)
		     nil #+nil(#_fftw_destroy_plan plan)))
		 (let* ((d (make-array (array-dimensions a)
				       :element-type '(complex double-float)))
			(d1 (make-array (* w h) :element-type (array-element-type d)
					:displaced-to d)))
		   (dotimes (i (length d1))
		     (setf (aref d1 i) (complex (aref c (* 2 i))       
						(aref c (+ 1 (* 2 i)))) ))
		   d))
	     (ccl:dispose-heap-ivector b)
	     (ccl:dispose-heap-ivector c)))))))

(defvar *bla* (make-array (list 512 512) :element-type '(unsigned-byte 16)))
(defvar *kbla* (make-array (list 512 512) :element-type '(unsigned-byte 16)))
(defvar *bla-sep* (make-array (list 128 128) :element-type '(complex double-float)))



#+nil
(defparameter *kbla* (scale-ub16 (ft *bla*)))

#+nil
(progn
  (setf *bla-sep* (ft (extract (ft *bla*)
			    :x (+ 256 -192) :y (+ 256 45) :w 128 :h 128)))
  nil)

(defun extract (a &key
                (x (floor (array-dimension a 1) 2))
                (y (floor (array-dimension a 0) 2)) 
                (w (next-power-of-two (min x y 
                                           (- (array-dimension a 1) x)
                                           (- (array-dimension a 0) y))))
                (h w))
  (let* ((b1 (make-array (* h w) :element-type (array-element-type a)
                         :initial-element 0))
         (b (make-array (list h w)
                        :element-type (array-element-type a)
                        :displaced-to b1))
         (ox (- x (floor w 2)))
         (oy (- y (floor h 2))))
    (assert (<= 0 ox))
    (assert (<= 0 oy))
    (assert (< (+ w ox) (array-dimension a 1)))
    (assert (< (+ h oy) (array-dimension a 0)))
    (dotimes (j h)
      (dotimes (i w)
        (setf (aref b j i)
              (aref a (+ j oy) (+ i ox)))))
    b))


(defun scale-ub16 (a)
  (declare (type (array (complex double-float) 2) a))
  (let* ((b (make-array (array-dimensions a) :element-type '(unsigned-byte 16)))
	 (b1 (make-array (reduce #'* (array-dimensions a)) :element-type '(unsigned-byte 16)
			 :displaced-to b))
	 (a1 (make-array (reduce #'* (array-dimensions a)) :element-type (array-element-type a)
			 :displaced-to a)))
    (dotimes (i (length b1))
      (setf (aref b1 i) (max 0 (min 65535 (floor (abs (aref a1 i)) 500)))))
    b))



;;ccl::*shared-libraries*
;;ccl::*eeps* ;; hash table with external functions

#+nil(progn (use-interface-dir :andor)
       (open-shared-library "libandor.so"))
#+nil
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
  (defun initialize (&optional (dir "/usr/local/etc/andor"))
    (with-cstrs ((cdir dir))
      (let ((r (#_Initialize cdir)))
	(unless (= #$DRV_SUCCESS r)
	  (break "Initialize returned ~a." r)))))
  (defun start-acquisition ()
    (assert (= #$DRV_SUCCESS (#_StartAcquisition))))
  (defun get-status ()
    (rlet ((status :int))
      (assert (= #$DRV_SUCCESS
		 (#_GetStatus
		  (%int-to-ptr (%address-of status)))))
      (%get-signed-long (%int-to-ptr (%address-of status)))))
  (defun set-read-mode (&optional (n 4))
    (assert (= #$DRV_SUCCESS (#_SetReadMode n))))
  (defun set-image (&key (hbin 1) (vbin 1) (hstart 1) (hend 1024) (vstart 1) (vend 1024))
    (assert (= #$DRV_SUCCESS (#_SetImage hbin vbin hstart hend vstart vend))))
  (defun set-exposure-time (&optional (time_s .1s0))
    (assert (= #$DRV_SUCCESS (#_SetExposureTime time_s))))
  (defun set-ad-channel (&optional (channel 0))
    (assert (= #$DRV_SUCCESS (#_SetADChannel channel))))
  (defun get-number-amp ()
    (rlet ((n :int))
      (assert (= #$DRV_SUCCESS
		 (#_GetNumberAmp
		  (%int-to-ptr (%address-of n)))))
      (%get-signed-long (%int-to-ptr (%address-of n)))))

  (defun get-number-available-images ()
    (rlet ((first :int)
	   (last :int))
      (let ((r (#_GetNumberAvailableImages
		(%int-to-ptr (%address-of first))
		(%int-to-ptr (%address-of last)))))
	(unless (= #$DRV_SUCCESS r)
	  (break "GetNumberAvailableImages returned ~a." r)))
      
      (list
       (%get-signed-long (%int-to-ptr (%address-of first)))
       (%get-signed-long (%int-to-ptr (%address-of last))))))


  (defun get-number-hs-speeds (&key (channel 1) (typ 1))
    (rlet ((speeds :int))
      (let ((r (#_GetNumberHSSpeeds channel typ 
				    (%int-to-ptr (%address-of
						  speeds)))))
	(unless (= #$DRV_SUCCESS r)
	  (break "GetNumberHSSpeeds returned ~a." r)))
      (%get-signed-long (%int-to-ptr (%address-of speeds)))))

  (defun get-hs-speed (&key (channel 1) (typ 1) (index 0))
    (rlet ((speeds :float))
      (let ((r (#_GetHSSpeed channel typ index
			     (%int-to-ptr (%address-of
					   speeds)))))
	(unless (= #$DRV_SUCCESS r)
	  (break "GetHSSpeed returned ~a." r)))
      (%get-single-float (%int-to-ptr (%address-of speeds)))))

  (defun set-hs-speed (&key (typ 1) (index 0))
    (let ((r (#_SetHSSpeed typ index)))
      (unless (= #$DRV_SUCCESS r)
	(break "SetHSSpeed returned ~a." r))))
  (defun abort-acquisition ()
    (assert (= #$DRV_SUCCESS (#_AbortAcquisition))))

  (defun get-acquisition-timings ()
    (rlet ((exposure :float)
	   (accumulate :float)
	   (kinetic :float))
      (let ((r (#_GetAcquisitionTimings (%int-to-ptr (%address-of exposure))
					(%int-to-ptr (%address-of accumulate))
					(%int-to-ptr (%address-of kinetic)))))
	(unless (= #$DRV_SUCCESS r)
	  (break "GetAcquisitionTimings returned ~a." r)))
      (list (%get-single-float (%int-to-ptr (%address-of exposure)))
	    (%get-single-float (%int-to-ptr (%address-of accumulate)))
	    (%get-single-float (%int-to-ptr (%address-of kinetic))))))

  (defun free-internal-memory ()
    (let ((r (#_FreeInternalMemory)))
      (unless (= #$DRV_SUCCESS r)
	(break "FreeInternalMemory returned ~a." r))))
  
  (defun get-number-ad-channels ()
    (rlet ((n :int))
      (assert (= #$DRV_SUCCESS
		 (#_GetNumberADChannels
		  (%int-to-ptr (%address-of n)))))
      (%get-signed-long (%int-to-ptr (%address-of n)))))

  (defun shut-down ()
    (assert (= #$DRV_SUCCESS (#_ShutDown))))
  (defun set-acquisition-mode (&optional (mode 1))
    (let ((r (#_SetAcquisitionMode mode)))
      (unless (= r #$DRV_SUCCESS)
	(break "SetAcquisitionMode returned ~a." r))))
  (defun get-acquired-data16 ()
    (let* ((w 1024)
	   (h 1024)
	   (n (* w h)))
      (multiple-value-bind (a ap) (make-heap-ivector n ;:unsigned-halfword 
						     '(unsigned-byte 16)
						     )
	(let ((r (#_GetAcquiredData16 ap n)))
	  (unless (= r #$DRV_SUCCESS)
	    (break "GetAcquiredData16 returned ~a." r)))
	(prog1 
	    (let* ((b1 (make-array n :element-type '(unsigned-byte 16)))
		   (b (make-array (list h w) :element-type '(unsigned-byte 16)
				  :displaced-to b1)))
	      (dotimes (i n)
		(setf (aref b1 i) (aref a i)))
	      b)
	  (dispose-heap-ivector a))))))




#+nil
(list (get-number-ad-channels)
      (get-number-amp)
      (get-number-hs-speeds :channel 1)
      (get-hs-speed :channel 1 :index 0))


#+nil
(set-hs-speed :typ 1 :index 0)

#+nil
(shut-down)



#+nil
(progn
  (initialize)
  (set-read-mode)
  (set-image :hend 1024 :vend 1024)
  (set-ad-channel 1)
  (set-exposure-time .00001)
  (set-acquisition-mode 1)
  (set-hs-speed :typ 1 :index 0)
  (get-acquisition-timings))

#+nil
(progn
 (defun get-time ()
   (rlet ((tv :timeval))
     (#_gettimeofday tv (%null-ptr))
     (+ (* 1000000 (pref tv :timeval.tv_sec))
	(pref tv :timeval.tv_usec))))
 (defun acquire-one-image ()
   (start-acquisition) 
   (unwind-protect 
	(progn
	  (let ((start (get-time)))
	    (loop while (and (/= #$DRV_IDLE (get-status))
			     (< (- (get-time) start) 1e6)) do
		 (sleep .01))
	    (format t "aquisition took ~a us~%" (- (get-time) start)))
	  (when (= #$DRV_IDLE (get-status))
	    (defparameter *bla* 
	      (get-acquired-data16))))
     (progn
       (when (= #$DRV_ACQUIRING (get-status))
	 (abort-acquisition))
       (free-internal-memory)))  
   (sleep .1)))


#+nil
(time
 (progn (acquire-one-image)
	(defparameter *kbla* (scale-ub16 (ft *bla*)))))

#+nil
(defparameter *kbla* (scale-ub16 (ft *bla*)))

(defparameter *fft-p* t)
#+nil
(ccl:process-run-function "fft" 
			  #'(lambda ()
			      (loop while *fft-p* do
				   (let ((f (ft *bla*)))
				     (defparameter *kbla* (scale-ub16 f))
				     (setf *bla-sep* 
					   (ft (extract f
							:x (+ 256 -192) 
							:y (+ 256 45)
							:w 128 :h 128)))
				     (dotimes (i 128)
				       (dotimes (j 128)
					 (setf (aref *bla-sep* j i)
					       (* (expt -1 (+ j i)) (aref *bla-sep* j i)))))))))
#+nil
(time
 (defparameter *kbla* (scale-ub16 (ft *bla*))))

(defparameter *acquire-p* t)
#+nil
(ccl:process-run-function "acquisition" 
			  #'(lambda ()
			      (loop while *acquire-p* do
				   (acquire-one-image)
				   )))

(defparameter *cap-num* 0)
#+nil
(ccl:process-run-function "acquisition" 
			  #'(lambda ()
			      (loop while *acquire-p* do
				   (capture-and-copy-frame)
				   (write-pgm (format nil "/media/home/martin/dat/bla~5,'0d.pgm" *cap-num*) *bla*)
				   (incf *cap-num*)
				   )))

#+nil
(capture-and-copy-frame)


(defvar *buffers* nil)
(progn
  (defun clear-frame ()
    (let* ((b1 (make-array (* 512 512) :element-type '(unsigned-byte 64)
			   :initial-element 0))
	   (b (make-array (list 512 512) :element-type '(unsigned-byte 64)
			 :displaced-to b1)))
      (defparameter *bla* b)))
  (defun capture-and-accumulate-frame ()
    (let* ((k (wait-and-read-frame))
	   (b1 (make-array (* 512 512) :element-type '(unsigned-byte 64)))
	   (b (make-array (list 512 512) :element-type '(unsigned-byte 64)
			 :displaced-to b1)))
      (declare (type (array (unsigned-byte 64) 2) b)
	       (type (array (unsigned-byte 64) 1) b1))
      (destructuring-bind  (ap length a a4) *buffers*
	(declare (type (array (unsigned-byte 8) 1) a)
		 (type (array (unsigned-byte 8) 4) a4))
	(destructuring-bind (kmax h w c) (array-dimensions a4)
	  (dotimes (i (min 512 w))
	   (dotimes (j (min 512 h))
	     (setf (aref b j i) (+ (aref b j i)
				   (aref a4 k (* 2 j) (* 2 i) 0)))))))
     (defparameter *bla* b))))

#+nil
(progn
  (write-mirror -70)
  (dotimes (i 20)
      (wait-and-read-frame))
  (clear-frame)
  (dotimes (i 10)
    (capture-and-accumulate-frame)))


#+nil
(clear-frame)

#+nil
(initialize)

#+nil
(get-status)

#+nil
(abort-acquisition)
#+nil
(set-read-mode)

#+nil
(set-image :hend 512 :vend 512)

#+nil
(set-exposure-time .01)

;; sudo mount --bind /dev/bus/ /proc/bus/      
;;cd /etc/udev/rules.d
;;sed -i -e s/SYSFS/ATTR/g *.rules
#+nil
(external-call "print_cam" :void)

#.(use-interface-dir :v4l2)

(progn
  (defvar *v4l-fd* nil)
  (defun v4l-open (&optional (fn "/dev/video1"))
    (unless *v4l-fd*
      (setf *v4l-fd* (with-cstrs ((fns fn))
		       (#_open fns (logior #$O_RDWR
					   #$O_NONBLOCK))))))
  (defun v4l-close ()
    (when *v4l-fd*
      (#_close *v4l-fd*))
    (setf *v4l-fd* nil))

  (defvar *buffers* nil)
  (defun v4l-allocate-buffers (&key (w 1600) (h 1200) (bytes-per-pixel 2))
    (unless *buffers*
      (setf *buffers* (let* ((number-buffers 4)
			     (length (* w h bytes-per-pixel))
			     (n (* number-buffers length))) 
			(multiple-value-bind (a ap) (make-heap-ivector n '(unsigned-byte 8))
			  (list ap length a (make-array (list number-buffers h w bytes-per-pixel)
							:element-type '(unsigned-byte 8)
							:displaced-to a)))))))
  (defun v4l-clear-buffers ()
    (when *buffers*
      (destructuring-bind (ap length a a4) *buffers*
	(dispose-heap-ivector a))
      (setf *buffers* nil)))
  (defun v4l-queue-buffers ()
    (destructuring-bind (ap length a a4) *buffers*
      (destructuring-bind (k h w c) (array-dimensions a4)
	(loop for i below k collect
	     (rletz ((buf :v4l2_buffer
			  :type #$V4L2_BUF_TYPE_VIDEO_CAPTURE
			  :memory #$V4L2_MEMORY_USERPTR
			  :index i
			  :m.userptr (+ (* i length) (%ptr-to-int ap))
			  :length length))
	       (let ((r (#_ioctl *v4l-fd* #$VIDIOC_QBUF :address buf)))
		 (assert (= 0 r))
		 r))))))
  (defun v4l-stream-on ()
    (multiple-value-bind (a ap) (make-heap-ivector 1 '(signed-byte 64))
      (setf (aref a 0) #$V4L2_BUF_TYPE_VIDEO_CAPTURE)
      (assert (= 0 (#_ioctl *v4l-fd* #$VIDIOC_STREAMON :address ap)))
      (dispose-heap-ivector a)))
  (defun v4l-stream-off ()
    (multiple-value-bind (a ap) (make-heap-ivector 1 '(signed-byte 64))
      (setf (aref a 0) #$V4L2_BUF_TYPE_VIDEO_CAPTURE)
      (assert (= 0 (#_ioctl *v4l-fd* #$VIDIOC_STREAMOFF :address ap)))
      (dispose-heap-ivector a)))
  (defun v4l-set-format (&key (w 1600) (h 1200))
    (rletz ((f :v4l2_format
	       :type #$V4L2_BUF_TYPE_VIDEO_CAPTURE
	       :fmt.pix.width w
	       :fmt.pix.height h
	       :fmt.pix.pixelformat #$V4L2_PIX_FMT_YUYV))
      (assert (= 0 (#_ioctl *v4l-fd* #$VIDIOC_S_FMT :address f)))
      (assert (= 0 (#_ioctl *v4l-fd* #$VIDIOC_G_FMT :address f)))
      (list (pref f :v4l2_format.fmt.pix.width)
	    (pref f :v4l2_format.fmt.pix.height)
	    :yuyv (= #$V4L2_PIX_FMT_YUYV
		     (pref f :v4l2_format.fmt.pix.pixelformat)))))
  (defun v4l-switch-to-user-pointers ()
    (rletz ((req :v4l2_requestbuffers ;; i want to read data with user pointers
		 :count 4
		 :type #$V4L2_BUF_TYPE_VIDEO_CAPTURE
		 :memory #$V4L2_MEMORY_USERPTR))
      (assert (= 0 (#_ioctl *v4l-fd* #$VIDIOC_REQBUFS :address req)))))
  (defun v4l-init (&key (fn "/dev/video1") (w 640) (h 480))
    (unless *v4l-fd*
      (v4l-open fn))
    (v4l-set-format :w w :h h)
    (v4l-switch-to-user-pointers)
    
    (v4l-clear-buffers)
    (v4l-allocate-buffers :w w :h h)
    (v4l-queue-buffers)
    (v4l-stream-on))
  (defun v4l-uninit ()
    (v4l-stream-off)
    (v4l-clear-buffers)
    (v4l-close))
  (defun read-frame ()
    (rletz ((buf :v4l2_buffer
		 :type #$V4L2_BUF_TYPE_VIDEO_CAPTURE
		 :memory #$V4L2_MEMORY_USERPTR))
      (assert (= 0 (#_ioctl *v4l-fd* #$VIDIOC_DQBUF :address buf)))
      (prog1 
	  (destructuring-bind (ap length a a4) *buffers*
	    (floor (- (pref buf :v4l2_buffer.m.userptr) (%ptr-to-int ap))
		   length))
	(assert (= 0 (#_ioctl *v4l-fd* #$VIDIOC_QBUF :address buf))))))
  (defun wait-and-read-frame ()
    (%stack-block ((fdset ccl::*fd-set-size*))
      (ccl::fd-zero fdset)
      (ccl::fd-set *v4l-fd* fdset)
      (rletz ((tv :timeval
		  :tv_sec 2
		  :tv_usec 0))
	(tagbody
	 :again
	 (let ((r (#_select (+ 1 *v4l-fd*) fdset 
			    (%null-ptr) (%null-ptr) tv)))
	   (ecase r
	     (1 (return-from wait-and-read-frame (read-frame)))
	     (0 (break "error: select timed out."))
	     (-1 (break "error: select returned ~a."
			(let ((err (ccl::%get-errno)))
			  (if (= err (- #$EINTR))
			      (go :again)
			      (list err (ccl::%strerror err)))))))))))))


#+nil
(v4l-uninit)
#+nil
(v4l-init :fn "/dev/video1" :w 800 :h 600)
#+nil
(unless *v4l-fd*
  (v4l-open "/dev/video0"))
#+nil
(v4l-set-format :w 640 :h 640)
#+nil
(progn (v4l-stream-off)
       (v4l-stream-on)
       (wait-and-read-frame))
#+nil
(progn
  (v4l-init :fn "/dev/video0")
  (loop for i below 10 collect (wait-and-read-frame)))

#+nil
(wait-and-read-frame)
#+nil
(capture-and-copy-frame)

(progn
  (defun capture-and-copy-frame ()
    (let* ((k (wait-and-read-frame))
	   (b1 (make-array (* 512 512) :element-type '(unsigned-byte 16)))
	   (b (make-array (list 512 512) :element-type '(unsigned-byte 16)
			  :displaced-to b1)))
      (declare (type (array (unsigned-byte 8) 2) b)
	       (type (array (unsigned-byte 8) 1) b1))
      (when k
	(destructuring-bind  (ap length a a4) *buffers*
	  (declare (type (array (unsigned-byte 8) 1) a)
		   (type (array (unsigned-byte 8) 4) a4))
	  (destructuring-bind (kmax h w c) (array-dimensions a4)
	    (dotimes (i (min 512 w))
	      (dotimes (j (min 512 h))
		(setf (aref b j i) (aref a4 k (* 1 j) (* 1 i) 0))))))
	(defparameter *bla* b))))


  #+nil
  (capture-and-copy-frame)

  (defun draw-texture (ap &key (w 512) (h 512))
    (let ((obj (first (gen-textures 1))))
      (bind-texture :texture-2d obj)
      (tex-parameter :texture-2d :texture-min-filter :linear)
      (tex-parameter :texture-2d :texture-mag-filter :linear)
      (tex-image-2d :texture-2d 0 :rgba w h 0
		    :green :unsigned-short ap)
      (enable :texture-2d)
      (with-primitive :quads
	(vertex 0 0) (tex-coord 0 0)
	(vertex 0 1) (tex-coord 0 1)
	(vertex 1 1) (tex-coord 1 1)
	(vertex 1 0) (tex-coord 1 0))
      (disable :texture-2d)
      (delete-textures (list obj))))
  
  (let ((rot 0))
    (defun draw ()
      (clear :color-buffer :depth-buffer)
      (color 1 1 1)
      (incf rot 1)
      (if (< 360 rot)
	  (setf rot 0))
      (with-pushed-matrix
	(translate -1.2 -2.4 -5)
	(rotate (* 3 rot) 1 0 0)
	(line-width 2)
	(let ((r .3))
	  (glut:wire-sphere r 15 7)
	  (color 0 0 0 1)
	  (glut:solid-sphere (* .99 r) 15 7)))

      (with-pushed-matrix
	(translate (sin (* 8 2 pi (/ rot 360))) 0 -2.2)
	(color 1 1 1)
	(rect .14 -3 -.14 3))
      (with-pushed-matrix 
	(translate -0.5 0 -1.4)
	(scale .5 .5 1)
	(color 1 0 0) (rect 0 0 1 .01)
	(color 0 1 0) (rect 0 0 .01 1))
      (let ((s .4))
       (with-pushed-matrix ;; reconstructed data amplitude
	 (translate -.5 0 -1.3)
	 (scale s s 1)
	 (color 1 1 1)
	 (let ((b1 (make-array (* 128 128) 
			       :element-type (array-element-type *bla-sep*)
			       :displaced-to *bla-sep*)))
	   (multiple-value-bind (a ap) (make-heap-ivector (* 128 128)
							  '(unsigned-byte 16))
	     (dotimes (i (* 128 128))
	       (setf (aref a i) (min 65535 (max 0 (floor (* .08 (abs (aref b1 i))))))))
	     (draw-texture ap :w 128 :h 128) 
	     (dispose-heap-ivector a))))
       (with-pushed-matrix ;; reconstructed data phase
	 (translate 0 0 -1.3)
	 (scale s s 1)
	 (color 1 1 1)
	 (let ((b1 (make-array (* 128 128) 
			       :element-type (array-element-type *bla-sep*)
			       :displaced-to *bla-sep*)))
	   (multiple-value-bind (a ap) (make-heap-ivector (* 128 128)
							  '(unsigned-byte 16))
	     (dotimes (i (* 128 128))
	       (setf (aref a i) (min 65535 (max 0 (floor
						   (* (/ 65535 (* 2 pi))
						      (+ pi (phase (aref b1 i)))))))))
	     (draw-texture ap :w 128 :h 128) 
	     (dispose-heap-ivector a)))))
      (with-pushed-matrix ;; fourier transform
	(translate -.5 0 -1.8)
	(scale 1 1 1)
	(color 1 1 1)
	(let ((b1 (make-array (* 512 512) :element-type '(unsigned-byte 16)
				:displaced-to *kbla*)))
	    (multiple-value-bind (a ap) (make-heap-ivector (* 512 512) '(unsigned-byte 16))
	      (dotimes (i (* 512 512))
		(setf (aref a i) (min 65535 (max 0 (* 256 100  (+ (aref b1 i)))))))
	      (draw-texture ap :w 512 :h 512) 
	      (dispose-heap-ivector a))))
      (with-pushed-matrix ;; camera image
	(translate -0.5 -1.0 -1.8)
	(scale 1 1 1)
	(color 1 1 1)
	(let ((b1 (make-array (* 512 512)
			      :element-type (array-element-type *bla*)
			      :displaced-to *bla*)))
	  (multiple-value-bind (a ap) (make-heap-ivector (* 512 512)
							 '(unsigned-byte 16))
	    (dotimes (i (* 512 512))
	      (setf (aref a i) (min 65535 (max 0 (* 256
						    (+ -0 (aref b1 i)))))))
	    (draw-texture ap :w 512 :h 512)
	    (dispose-heap-ivector a))))
      (with-pushed-matrix  ;; histogram
	(translate -.47 -.07 -1.3)
	(scale (/ 512.0) .005 1)
	(color 1 1 1)
	(let ((hist (make-array 256 :element-type 'fixnum)))
	  (destructuring-bind (h w) (array-dimensions *bla*)
	    (dotimes (i w)
	      (dotimes (j h)
		(incf (aref hist (min 255 (aref *bla* j i)))))))
	  (with-primitive :lines 
	    (vertex 0 10)
	    (vertex 256 10)
	    (dotimes (i 256)
	      (vertex i 0)
	      (vertex i (log (+ 1 (aref hist i))))))))
      (flush)
      (sleep (/ 61.3))
					;   (cl-opengl-bindings:wait-sync )
      (glut:post-redisplay))))


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
