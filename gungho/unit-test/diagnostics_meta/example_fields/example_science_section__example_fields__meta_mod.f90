!-------------------------------------------------------------------------------
! (C) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!
!               Sample Module Diagnostic Meta File
!
! This is a sample file that illustrates how you define meta data for a field
! THIS IS TO BE UPDATED WITH DETAILS ABOUT WHERE THESE FILES ARE TO BE LOCATED,
! ONCE THIS IS DECIDED.
!
! (*/source/diagnostics_meta/*_meta_mod.f90)
!
!-------------------------------------------------------------------------------
module example_science_section__example_fields__meta_mod

  use meta_data_mod,                  only: meta_data_type
  use constants_mod,                  only: real_type, r_def, i_def
  !> Only import the dimensions that you will actually be using
  use vertical_dimensions_mod,        only: model_height_dimension, &
                                            model_depth_dimension, &
                                            fixed_height_dimension
  use misc_meta_data_mod,             only: misc_meta_data_type
  !> Only import the function spaces that you will actually be using
  use fs_continuity_mod,              only: W2, W3, Wtheta
  !> Only import the time steps that you will actually be using
  use time_step_enum_mod,             only: STANDARD_TIMESTEP
  !> Only import the interpolation methods that you will actually be using
  use interpolation_enum_mod,         only: BILINEAR
  !> Only import the levels that you will actually be using
  use levels_enum_mod,                only: TOP_WET_LEVEL

  implicit none

  private

  type, public :: example_science_section__example_fields__meta_type

    !> Declare the name of your fields here
    type(meta_data_type), public :: &
      eastward_wind, &
      rate_of_increase_rain_mass_due_to_autoconv_from_liquid_cloud, &
      air_potential_temperature, &
      lwe_thickness_of_moisture_content_of_soil_layer, &
      surface_altitude, &
      air_temperature_over_tiles, &
      low_type_cloud_area_fraction

    end type example_science_section__example_fields__meta_type

  interface example_science_section__example_fields__meta_type
    module procedure example_science_section__example_fields__meta_constructor
  end interface

contains

  !>@brief Creates meta_data_type objects for a specific section of science
  function example_science_section__example_fields__meta_constructor() result(self)
    implicit none

    type(example_science_section__example_fields__meta_type) :: self

    !> Example field using a height vertical dimension using model levels
    !> If no arguments are present, it will default to top atmospheric level
    !> and bottom atmospheric level. This field also uses a standard name,
    !> which is optional
    self%eastward_wind = meta_data_type(&
      unique_id = "example_fields__eastward_wind", &
      units = "m s-1", &
      function_space = W2, &
      order = 0, &
      io_driver = "", &
      trigger = "_checksum: true;", &
      description = "u component of wind on u pts on native c grid.", &
      data_type = REAL_TYPE, &
      time_step = STANDARD_TIMESTEP, &
      recommended_interpolation = BILINEAR, &
      packing = 0, &
      vertical_dimension = model_height_dimension(), &
      standard_name = "eastward_wind", &
      misc_meta_data = [misc_meta_data_type("test_name_1","test_value_1"), &
                        misc_meta_data_type("test_name_2","test_value_2")])


    !> Example of a field using a model height dimension and supplying one
    !> of its arguments
    self%rate_of_increase_rain_mass_due_to_autoconv_from_liquid_cloud = meta_data_type(&
      unique_id = "example_fields__rate_of_increase_of_rain_mass_due_to_"// &
                  "autoconversion_from_liquid_cloud", &
      units = "kg kg-1 s-1", &
      function_space = Wtheta, &
      order = 0, &
      io_driver = "", &
      trigger = "_checksum: true;", &
      description = "This is a microphysical process transfer rate. "// &
                    "Outputting all microphysics", &
      data_type = REAL_TYPE, &
      time_step = STANDARD_TIMESTEP, &
      recommended_interpolation = BILINEAR, &
      packing = 0, &
      vertical_dimension = model_height_dimension(top = TOP_WET_LEVEL))

    !> This field uses a different function space
    self%air_potential_temperature = meta_data_type(&
      unique_id = "example_fields__air_potential_temperature", &
      units = "K", &
      function_space = Wtheta, &
      order = 0, &
      io_driver = "", &
      trigger = "_checksum: true;", &
      description = "Potential temperature on p points on native c grid.", &
      data_type = REAL_TYPE, &
      time_step = STANDARD_TIMESTEP, &
      recommended_interpolation = BILINEAR, &
      packing = 0, &
      vertical_dimension = model_height_dimension(), &
      standard_name = "air_potential_temperature")

    !> Example field using a depth vertical dimension using model levels
    !> If no arguments are present, it will default to top soil level and
    !> bottom soil level
    self%lwe_thickness_of_moisture_content_of_soil_layer = meta_data_type(&
      unique_id = "example_fields__lwe_thickness_of_moisture_content_of_"// &
                  "soil_layer", &
      units = "kg m-2", &
      function_space = W3, &
      order = 0, &
      io_driver = "", &
      trigger = "_checksum: true;", &
      description = "Total (frozen+unfrozen) soil moisture content in a "// &
                    "soil layer (kg/m2). This example field will be dealt"// &
                    " with using negative height values (Phase 3)", &
      data_type = REAL_TYPE, &
      time_step = STANDARD_TIMESTEP, &
      recommended_interpolation = BILINEAR, &
      packing = 0, &
      vertical_dimension = model_depth_dimension(), &
      standard_name = "lwe_thickness_of_moisture_content_of_soil_layer")

    !> Example field using a height vertical dimension on fixed levels
    !> The level definition should be passed as an array of
    !> floating point numbers
    self%surface_altitude = meta_data_type(&
      unique_id = "example_fields__surface_altitude", &
      units = "m", &
      function_space = Wtheta, &
      order = 0, &
      io_driver = "", &
      trigger = "_checksum: true;", &
      description = 'The surface called "surface" means the lower boundary '// &
                    'of the atmosphere. Altitude is the (geometric) height '// &
                    'above the geoid, which is the reference geopotential '// &
                    'surface. The geoid is similar to mean sea level.', &
      data_type = REAL_TYPE, &
      time_step = STANDARD_TIMESTEP, &
      recommended_interpolation = BILINEAR, &
      packing = 0, &
      vertical_dimension = fixed_height_dimension(level_definition = REAL([0.0], r_def)), &
      standard_name = "surface_altitude")

    !> Another example of a field using a height vertical
    !> dimension on fixed levels.
    !> This field will require an additional non-spatial dimension to
    !> define the tiles. This is yet to be implemented
    self%air_temperature_over_tiles = meta_data_type(&
      unique_id = "example_fields__air_temperature_over_tiles", &
      units = "K", &
      function_space = WTheta, &
      order = 0, &
      io_driver = "", &
      trigger = "_checksum: true;", &
      description = "1.5M TEMPERATURE OVER TILES", &
      data_type = REAL_TYPE, &
      time_step = STANDARD_TIMESTEP, &
      recommended_interpolation = BILINEAR, &
      packing = 0, &
      vertical_dimension = fixed_height_dimension(level_definition = REAL([1.5], r_def)), &
      standard_name = "air_temperature")

    !> Example of a field using a fixed height dimension
    self%low_type_cloud_area_fraction = meta_data_type(&
      unique_id = "example_fields__low_type_cloud_area_fraction", &
      units = "1", &
      function_space = W3, &
      order = 0, &
      io_driver = "", &
      trigger = "_checksum: true;", &
      description = "Low cloud amount. Phase 3", &
      data_type = REAL_TYPE, &
      time_step = STANDARD_TIMESTEP, &
      recommended_interpolation = BILINEAR, &
      packing = 0, &
      vertical_dimension = fixed_height_dimension(level_definition = REAL([111.0, 1949.0], r_def)), &
      standard_name = "low_type_cloud_area_fraction")

  end function example_science_section__example_fields__meta_constructor
end module example_science_section__example_fields__meta_mod