#.(progn
    (use-interface-dir :fftw3)
    (open-shared-library "libfftw3.so"))

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
      (dotimes (j h)
        (dotimes (i w)
         (write-byte (aref a j i) s))))))
(defun read-pgm (filename)
  (declare ((or pathname string) filename)
	   #+sbcl (values (or (simple-array (unsigned-byte 8) 2)
			      (simple-array (unsigned-byte 16) 2)) &optional))
  (with-open-file (s filename)
    (unless (equal (symbol-name (read s)) "P5")
      (error "no PGM file"))
    (let* ((w (read s))
         (h (read s))
         (grays (read s))
         (pos (file-position s)))
      (declare ((integer 0 65535) grays w h))
      (let* ((type (if (<= grays 255)
                 '(unsigned-byte 8)
                 '(unsigned-byte 16)))
         (data (make-array (list h w)
                         :element-type type))
         (data-1d (make-array (* h w)
                                 :element-type type
                                 :displaced-to data)))
        (with-open-file (s2 filename :element-type type)
         (file-position s2 pos)
         (read-sequence data-1d s2))
               data))))

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


(defun scale-ub8 (a &key (fun #'abs) (scale .01d0) (offset 0d0))
  (declare (type (array (complex double-float) 2) a))
  (let* ((b (make-array (array-dimensions a) :element-type '(unsigned-byte 8)))
	 (b1 (make-array (reduce #'* (array-dimensions a)) :element-type '(unsigned-byte 8)
			 :displaced-to b))
	 (a1 (make-array (reduce #'* (array-dimensions a)) :element-type (array-element-type a)
			 :displaced-to a)))
    (dotimes (i (length b1))
      (setf (aref b1 i) (max 0 (min 255 (floor (* scale (- (funcall fun (aref a1 i)) offset)))))))
    b))

(array-dimension (make-array (list 3 2)) 1)


(defun cat (&rest r)
  (let* ((h (loop for e in r maximize (array-dimension e 0)))
	 (w (loop for e in r sum (array-dimension e 1)))
	 (c (make-array (list h w) :element-type '(unsigned-byte 8)))
	 (w0 0))
    (loop for e in r do
	 (destructuring-bind (h w) (array-dimensions e)
	   (dotimes (i w)
	     (dotimes (j h)
	       (setf (aref c j (+ w0 i)) (aref e j i))))
	   (incf w0 w)))
    c))

#+nil
(defparameter *bla*
  (write-pgm "/dev/shm/q2/o.pgm"
   (cat (read-pgm (format nil "/media/home/martin/dat/bla~5,'0d.pgm" 1))
	(read-pgm (format nil "/media/home/martin/dat/bla~5,'0d.pgm" 2)))))

(defparameter *bla-sep-old* (make-array (list 128 128) :element-type 'double-float))
(defparameter *bla-sep* (make-array (list 128 128) :element-type 'double-float))
(defparameter *dphi-sum* nil)

#+nil
(progn
  (with-open-file (p "/dev/shm/q2/dphi-sum.gp" :direction :output
		    :if-exists :supersede :if-does-not-exist :create)
    (format p "plot \"dphi-sum.dat\" u 1:2 w l; pause -1;"))
 (with-open-file (p "/dev/shm/q2/dphi-sum.dat" :direction :output
		    :if-exists :supersede :if-does-not-exist :create)
   (loop for (i f) in (reverse *dphi-sum*) do
	(format p "~d ~g~%" i f))))

#+nil
(with-open-file (log "/dev/shm/q2/log" :direction :output
		     :if-exists :supersede :if-does-not-exist :create)
 (loop for i below 10812 do
      (let* ((num (format nil "~5,'0d" i))
	     (img (read-pgm (format nil "/media/home/martin/dat/bla~a.pgm" num)))
	     (f (ft img)))
	(format log "~a~%" i) (force-output log)
	(defparameter *bla-sep* 
	  (ft (extract f
		       :x (+ 256 45)
		       :y (+ 256 -192) 
		       :w 128 :h 128)))
	(dotimes (i 128)
	  (dotimes (j 128)
	    (setf (aref *bla-sep* j i)
		  (* (expt -1 (+ j i)) (aref *bla-sep* j i)))))
	(write-pgm (format nil "/dev/shm/q2/o~a.pgm" num)
		   (cat (extract img :w 128 :h 128) 
			(scale-ub8 (extract f
					    :x (+ 256 45)
					    :y (+ 256 -192) 
					    :w 128 :h 128) :scale .003)
			(scale-ub8 *bla-sep* :scale .00003)
			(scale-ub8 *bla-sep* :fun #'phase :offset (- pi) :scale (* 255 .5 (/ pi)))
			(when *bla-sep-old*
			  (let ((dphi (make-array (list 128 128) :element-type 'double-float)))
			    (dotimes (i 128)
			      (dotimes (j 128)
				(setf (aref dphi j i) (if (< (abs (aref *bla-sep* j i)) 4d4)
							  0d0
							  (- (imagpart (/ (- (aref *bla-sep* j i)
									     (aref *bla-sep-old* j i))
									  (aref *bla-sep* j i))))))))
			    (let ((sum 0d0))
			      (dotimes (i 128)
				(dotimes (j 128)
				  (incf sum (* (aref dphi j i) (aref dphi j i)))))
			      (push (list i sum) *dphi-sum*))
			    (scale-ub8 dphi :fun #'identity :offset -.9 :scale 100.0)))))
	(defparameter *bla-sep-old* *bla-sep*))))


