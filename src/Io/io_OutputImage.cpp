// See LICENSE_CELLO file for license and copyright information

/// @file     io_OutputImage.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Mar 17 11:14:18 PDT 2011
/// @brief    Implementation of the OutputImage class

#include "cello.hpp"

#include "io.hpp"

//----------------------------------------------------------------------

OutputImage::OutputImage() throw ()
  : Output(),
    image_(0),
    image_size_x_(0),
    image_size_y_(0),
    png_(0)

{
    map_r_.resize(2);
    map_g_.resize(2);
    map_b_.resize(2);
    map_r_[0] = 0.0;
    map_g_[0] = 0.0;
    map_b_[0] = 0.0;
    map_r_[1] = 1.0;
    map_g_[1] = 1.0;
    map_b_[1] = 1.0;

    // Only root process writes
#ifdef CONFIG_USE_CHARM
    process_stride_ = CkNumPes();
#endif
}

//----------------------------------------------------------------------

OutputImage::~OutputImage() throw ()
{
}

//======================================================================

#ifdef CONFIG_USE_CHARM

void OutputImage::open (const Hierarchy * hierarchy, int cycle, double time) throw()
{
  // create process image and clear it

  int nxm,nym,nzm;
  hierarchy->patch(0)->size (&nxm, &nym, &nzm);

  image_create_(nxm,nym);

  int ip = CkMyPe();

  if (is_writer_(ip)) create_("hierarchy",nxm,nym,cycle,time);

}

//----------------------------------------------------------------------

void OutputImage::close () throw()
{
  close_();
}

//----------------------------------------------------------------------

void OutputImage::block (const Block * block) throw()
{
  //  write(block,
  
  //   incorporate block data into process data
}

#endif

//----------------------------------------------------------------------

void OutputImage::write
(
 const FieldDescr * field_descr,
 Hierarchy * hierarchy,
 int cycle,
 double time,
 bool root_call
  ) throw()
{

  if (root_call) {

    int nxm,nym,nzm;
    hierarchy->patch(0)->size (&nxm, &nym, &nzm);

    create_("hierarchy",nxm, nym, cycle,time);
  }

  ItPatch it_patch (hierarchy);
  while (Patch * patch = ++it_patch) {
    // NO OFFSET: ASSUMES ROOT PATCH
    write (field_descr, patch, cycle,time, false,  0,0,0);
  }

  if (root_call) close_();

}

//----------------------------------------------------------------------

void OutputImage::write
(
 const FieldDescr * field_descr,
 Patch * patch,
 int cycle,
 double time,
 bool root_call,
 int ix0,
 int iy0,
 int iz0
 ) throw()
{
  if (root_call) {
    int nxp, nyp;
    patch->size (&nxp, &nyp);
    create_("patch",nxp,nyp,cycle,time);
  }

#ifdef CONFIG_USE_CHARM

  if (patch->blocks_allocated()) {
    CkPrintf ("%s:%d Output blocks in Patch %p\n",__FILE__,__LINE__,patch);
    //    patch->blocks()
  }
    
#else

  ItBlock it_block (patch);
  Block * block;

  while ((block = ++it_block)) {

    FieldBlock * field_block = block->field_block();

    int nxb,nyb,nzb;
    field_block->size(&nxb,&nyb,&nzb);

    int ix,iy,iz;
    block->index_patch(&ix,&iy,&iz);

    write (field_descr, block, cycle,time,false,
	   ix0+ix*nxb,
	   iy0+iy*nyb,
	   iz0+iz*nzb);
  }

#endif

  if (root_call) close_();

}

//----------------------------------------------------------------------

void OutputImage::write
(
 const FieldDescr * field_descr,
 Block * block,
 int cycle,
 double time,
 bool root_call,
 int ix0,
 int iy0,
 int iz0
) throw()
{
  FieldBlock * field_block = block->field_block();

  // Get block size
  int nx,ny,nz;
  field_block->size(&nx,&ny,&nz);

  // Index of (only) field to write
  it_field_->first();
  int index = it_field_->value();

    // Get ghost depth

  int gx,gy,gz;
  field_descr->ghosts(index,&gx,&gy,&gz);

  // Get array dimensions
  int ndx,ndy,ndz;
  ndx=nx+2*gx;
  ndy=ny+2*gy;
  ndz=nz+2*gz;

  if (root_call) create_("block",nx,ny,cycle,time);

  // Add block contribution to image

  char * field_unknowns = field_block->field_unknowns(field_descr,index);
    
  if (field_descr->precision(index) == precision_single) {
    image_reduce_ (((float *) field_unknowns),
		   ndx,ndy,ndz, 
		   nx,ny,nz,
		   ix0,iy0,iz0,
		   axis_z, reduce_sum);
  } else {
    image_reduce_ (((double *) field_unknowns),
		   ndx,ndy,ndz, 
		   nx,ny,nz,
		   ix0,iy0,iz0,
		   axis_z, reduce_sum);
  }

  if (root_call) close_();

}

//======================================================================

void OutputImage::image_set_map
(int n, double * map_r, double * map_g, double * map_b) throw()
{
  map_r_.resize(n);
  map_g_.resize(n);
  map_b_.resize(n);

  for (int i=0; i<n; i++) {
    map_r_[i] = map_r[i];
    map_g_[i] = map_g[i];
    map_b_[i] = map_b[i];
  }
}

//----------------------------------------------------------------------

void OutputImage::png_create_ (std::string filename, int mx, int my) throw()
{
  png_ = new pngwriter(mx,my,0,filename.c_str());
}

//----------------------------------------------------------------------

void OutputImage::png_close_ () throw()
{
  png_->close();
  delete png_;
  png_ = 0;
}

//----------------------------------------------------------------------

void OutputImage::image_create_ (int mx, int my) throw()
{
  image_size_x_ = mx;
  image_size_y_ = my;

  image_ = new double [mx*my];

  for (int i=0; i<mx*my; i++) image_[i] = 0.0;
}

//----------------------------------------------------------------------

void OutputImage::image_close_ (double min, double max) throw()
{

  // simplified variable names

  int mx = image_size_x_;
  int my = image_size_y_;
  int m  = mx*my;

  // Adjust min and max bounds if needed

   for (int i=0; i<m; i++) {
     min = MIN(min,image_[i]);
     max = MAX(max,image_[i]);
   }

  // loop over pixels (ix,iy)

  for (int ix = 0; ix<mx; ix++) {

    for (int iy = 0; iy<my; iy++) {

      int i = ix + mx*iy;

      double value = image_[i];

      double r = 1.0, g = 0, b = 0;

      if (min <= image_[i] && image_[i] <= max) {

	// map v to lower colormap index
	size_t k = (map_r_.size() - 1)*(value - min) / (max-min);

	// prevent k == map_.size()-1, which happens if value == max
	if (k > map_r_.size() - 2) k = map_r_.size()-2;

	// linear interpolate colormap values
	double lo = min +  k   *(max-min)/(map_r_.size()-1);
	double hi = min + (k+1)*(max-min)/(map_r_.size()-1);

	// should be in bounds, but force if not due to rounding error
	if (value < lo) value = lo;
	if (value > hi) value = hi;

	// interpolate colormap

	double ratio = (value - lo) / (hi-lo);

	r = (1-ratio)*map_r_[k] + ratio*map_r_[k+1];
	g = (1-ratio)*map_g_[k] + ratio*map_g_[k+1];
	b = (1-ratio)*map_b_[k] + ratio*map_b_[k+1];
      }

      // Plot pixel, red if out of range
      png_->plot(ix+1, iy+1, r,g,b);
    }
  }      

  delete [] image_;
  image_ = 0;
}

//----------------------------------------------------------------------

void OutputImage::create_ 
(
 std::string write_level,
 int nx, int ny,
 int cycle,
 double time
 ) throw()
{
  // Open file if writing a single block
  std::string file_name = expand_file_name(cycle,time);

  // Create image 
  image_create_(nx,ny);

  // Create png object
  png_create_(file_name,nx,ny);

  Monitor::instance()->print ("[Output] writing %s image file %s", 
			      write_level.c_str(),file_name.c_str());
  
}

//----------------------------------------------------------------------

void OutputImage::close_ () throw()
{
  double min=0.0;
  double max=1.0;
  image_close_(min,max);
  png_close_();
}
