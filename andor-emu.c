// gcc -shared andor-emu.c -o libandor-emu.so 
#include "/usr/local/include/atmcdLXd.h"
#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>

enum { HANDLE = 100 };

typedef struct Camera Camera;
struct Camera {
  unsigned int width,height,bpp;//offx,offy,binx,biny;
  float exptime;
  unsigned short *buf;
  int initialized_p, acquiring_p, data_new_p;
  struct timeval tv;
};

struct Camera cam;

void print_cam()
{
  printf("w=%d, h=%d, bpp=%d, exptime=%g, buf=0x%lx\n",
	 cam.width,cam.height,cam.bpp,cam.exptime,(unsigned long)cam.buf);
  printf("initialized=%d, acquiring=%d, data_new=%d\n",
	 cam.initialized_p,cam.acquiring_p,cam.data_new_p);
  printf("tv=%ld %ld\n",
	 cam.tv.tv_sec,cam.tv.tv_usec);
  fflush(0);
}

unsigned int GetAvailableCameras(at_32 * totalCameras)
{
  *totalCameras = 1;
  return DRV_SUCCESS;
}

unsigned int GetCameraHandle(at_32 cameraIndex, at_32 * cameraHandle)
{
  if(cameraIndex == 0){
    *cameraHandle=HANDLE;
    return DRV_SUCCESS;
  }
  return DRV_P1INVALID;
}

unsigned int SetCurrentCamera(at_32 camera_handle)
{
  if(camera_handle == HANDLE)
    return DRV_SUCCESS;
  return DRV_P1INVALID;
}

unsigned int Initialize(char*s)
{
  if(s==0){}
  if(!cam.initialized_p){
    cam.width = 512;
    cam.height = 512;
    cam.bpp = 16;
    cam.exptime = .1;
    cam.buf = malloc(cam.width*cam.height*cam.bpp/8);
    cam.initialized_p = 1;
    cam.acquiring_p=0;
    cam.data_new_p=0;
  }
  return DRV_SUCCESS;
}

unsigned int StartAcquisition()
{
  gettimeofday(&(cam.tv),0);
  cam.acquiring_p = 1;
  cam.data_new_p = 0;
  return DRV_SUCCESS;
}

unsigned int ShutDown()
{
  cam.initialized_p=0;
  free(cam.buf);
  return DRV_SUCCESS;
}

long long get_time_since_acquisition_start()
{
      struct timeval tv,diff;
      gettimeofday(&tv,0);
      diff.tv_sec = tv.tv_sec - cam.tv.tv_sec;
      diff.tv_usec = tv.tv_usec - cam.tv.tv_usec;
      return 1000000*diff.tv_sec+diff.tv_usec;
}

unsigned int GetStatus(int*status)
{
  if(cam.initialized_p){
    if(cam.acquiring_p){
      if(get_time_since_acquisition_start()<cam.exptime){
	*status=DRV_ACQUIRING;
      } else {
	unsigned int i;
	for(i=0;i<cam.width*cam.height;i++)
	  cam.buf[i] = (unsigned short)(100.0+30.0*drand48());
	cam.data_new_p = 1;
	cam.acquiring_p = 0;
	*status = DRV_IDLE;
      }
      return DRV_SUCCESS;
    }
    *status = DRV_IDLE;
    return DRV_SUCCESS;
  }
  return DRV_NOT_INITIALIZED;
}

unsigned int GetAcquiredData16(unsigned short*a, at_u32 n)
{
  // a must have been allocated by the user
  if(a == 0)
    return DRV_P1INVALID;
  if(n != cam.width*cam.height)
    return DRV_P2INVALID;
  int status;
  GetStatus(&status);
  switch(status){
  case DRV_ACQUIRING: return DRV_ACQUIRING;
  case DRV_NOT_INITIALIZED: return DRV_NOT_INITIALIZED;
  }
  if(!cam.data_new_p)
    return DRV_NO_NEW_DATA;
  int i;
  for(i=0;i<n;i++)
    a[i]=cam.buf[i];
  cam.data_new_p=0;
  return DRV_SUCCESS;
}

/* more things i should implement
circular buffer and continuous acquisition mode
SetAcquisitionMode
SetExposureTime
SetTriggerMode
GetAcquisitionTimings
*/
