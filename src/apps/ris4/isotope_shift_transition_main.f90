! Combine two RIS4 state reports and nuclear data into a line isotope shift.
program isotope_shift_transition_main
    ! Module dependencies
    use calculation_state_files, only: read_state_pair, isotope_shift_input_path, lsj_label_input_path
    use, intrinsic :: iso_fortran_env, only: int64, real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use run_timer, only: begin_run_timer, finish_run_timer
    implicit none (type, external)
    integer :: isotope_shift_output_unit
    integer(int64) :: start_clock_count
    integer :: i,j,ii,jj,lsjpresent,lmax,lmax1,lmax2
    integer :: first_state_level_count,second_state_level_count,nveclsj1,nveclsj2,level1(500),level2(500),levellsj1(500), &
        & levellsj2(500)
    integer :: relflag,efacflag

    real(real64) :: energy1(500),nms1(500),sms1(500),ms1(500),fs1(6,500),fshift1(500),fshift11(500)
    real(real64) :: trms_1,t_1,w_1,b20_1,b40_1, radmom1(4)
    real(real64) :: trms_2,t_2,w_2,b20_2,b40_2, radmom2(4)

    real(real64) :: energy2(500),nms2(500),sms2(500),ms2(500),fs2(6,500),fshift2(500),fshift22(500)
    real(real64) :: tenergy, nmslaw, nmsconst, m1, m2, dr2, masshift,fieldshift,fieldshift22,fieldshift1(500),fieldshift2(500)
    real(real64) :: isotopeshift
    real(real64) :: f0, f2, f4, f6, f0res, f2res, linfac
    real(real64) :: fieldshiftres, fieldshift22res, linfac2

    real(real64), parameter :: hartree_to_wavenumber_cm = 219474.63136320_real64   ! CODATA 2018
    real(real64), parameter :: mev2ghz = 241.798926_real64
    real(real64) :: inverse_mass_difference, radius_difference
    real(real64) :: z

    character(len=101) :: state1, state2
    character(:), allocatable :: file1, file2, file3, file5, file6
    character(len=11) :: spinpar1(500),spinpar2(500)
    character(len=4) :: typeflag,spinlsj1(500),spinlsj2(500)
    character(len=3) :: unit
    character(len=1) :: def,ci,rel,efac,soph,parlsj1(500),parlsj2(500)
    character(len=200) :: csflabel1(500),csflabel2(500),csflabel1scal,csflabel2scal
    character(len=200) :: headermpghz,headermpmhz,headermpmev
    character(len=200) :: headermpghzlsj,headermpmhzlsj,headermpmevlsj
    character(len=200) :: headerfpghz,headerfpmhz,headerfpmev
    character(len=200) :: headerfpghzlsj,headerfpmhzlsj,headerfpmevlsj
    character(len=200) :: headerfp2ghz,headerfp2mhz,headerfp2mev
    character(len=200) :: headerfp2ghzlsj,headerfp2mhzlsj,headerfp2mevlsj
    character(len=200) :: headereghz,headeremhz,headeremev
    character(len=200) :: headereghzlsj,headeremhzlsj,headeremevlsj
    character(len=50) :: format1a,format1b,format1b2,format1c,format1c2,format2
    character(len=50) :: format1alsj,format1blsj,format1b2lsj,format1clsj,format1c2lsj,format2lsj
    character(len=200) :: uplevel,lowlevel

    logical :: exist1,exist2

    call begin_run_timer (start_clock_count, "FICAL")
    headermpghz = '(1X,"Upper level",5X,"Lower level",5X,"Energy (cm-1)",3X,"NMS-S (GHz u)" &
        &,6X,"NMS (GHz u)",6X,"SMS (GHz u)",8X,"MS (GHz u)")'
    headermpghzlsj = '(1X,a,5X,a,5X,"Energy (cm-1)",3X,"NMS-S (GHz u)" &
        &,6X,"NMS (GHz u)",6X,"SMS (GHz u)",8X,"MS (GHz u)")'
    headermpmhz = '(1X,"Upper level",5X,"Lower level",5X,"Energy (cm-1)",3X,"NMS-S (MHz u)" &
        &,6X,"NMS (MHz u)",6X,"SMS (MHz u)",8X,"MS (MHz u)")'
    headermpmhzlsj = '(1X,a,5X,a,5X,"Energy (cm-1)",3X,"NMS-S (MHz u)" &
        &,6X,"NMS (MHz u)",6X,"SMS (MHz u)",8X,"MS (MHz u)")'
    headermpmev = '(1X,"Upper level",5X,"Lower level",5X,"Energy (cm-1)",3X,"NMS-S (meV u)" &
        &,6X,"NMS (meV u)",6X,"SMS (meV u)",8X,"MS (meV u)")'
    headermpmevlsj = '(1X,a,5X,a,5X,"Energy (cm-1)",3X,"NMS-S (meV u)" &
        &,6X,"NMS (meV u)",6X,"SMS (meV u)",8X,"MS (meV u)")'

    headerfpghz = '(1X,"Upper level",5X,"Lower level",5X,"Energy (cm-1)" &
        &,3X,"F0 (GHz fm-2)",5X,"F2 (GHz fm-4)",5X,"F4 (GHz fm-6)",5X,"F6 (GHz fm-8)",5X,"Lin. dep. (per mille)")'
    headerfpghzlsj = '(1X,a,5X,a,5X,"Energy (cm-1)" &
        &,3X,"F0 (GHz fm-2)",5X,"F2 (GHz fm-4)",5X,"F4 (GHz fm-6)",5X,"F6 (GHz fm-8)",5X,"Lin. dep. (per mille)")'
    headerfpmhz = '(1X,"Upper level",5X,"Lower level",5X,"Energy (cm-1)" &
        &,3X,"F0 (MHz fm-2)",5X,"F2 (MHz fm-4)",5X,"F4 (MHz fm-6)",5X,"F6 (MHz fm-8)")'
    headerfpmhzlsj = '(1X,a,5X,a,5X,"Energy (cm-1)" &
        &,3X,"F0 (MHz fm-2)",5X,"F2 (MHz fm-4)",5X,"F4 (MHz fm-6)",5X,"F6 (MHz fm-8)")'
    headerfpmev = '(1X,"Upper level",5X,"Lower level",5X,"Energy (cm-1)" &
        &,3X,"F0 (meV fm-2)",5X,"F2 (meV fm-4)",5X,"F4 (meV fm-6)",5X,"F6 (meV fm-8)")'
    headerfpmevlsj = '(1X,a,5X,a,5X,"Energy (cm-1)" &
        &,3X,"F0 (meV fm-2)",5X,"F2 (meV fm-4)",5X,"F4 (meV fm-6)",5X,"F6 (meV fm-8)")'

    headerfp2ghz = '(1X,"Upper level",5X,"Lower level",5X,"Energy (cm-1)" &
        &,1X,"F0VED0 (GHz fm-2)",1X,"F0VED1 (GHz fm-4)",3X,"Lin. dep. (per mille)")'
    headerfp2ghzlsj = '(1X,a,5X,a,5X,"Energy (cm-1)" &
        &,1X,"F0VED0 (GHz fm-2)",1X,"F0VED1 (GHz fm-4)",3X,"Lin. dep. (per mille)")'
    headerfp2mhz = '(1X,"Upper level",5X,"Lower level",5X,"Energy (cm-1)" &
        &,1X,"F0VED0 (MHz fm-2)",1X,"F0VED1 (MHz fm-4)")'
    headerfp2mhzlsj = '(1X,a,5X,a,5X,"Energy (cm-1)" &
        &,1X,"F0VED0 (MHz fm-2)",1X,"F0VED1 (MHz fm-4)")'
    headerfp2mev = '(1X,"Upper level",5X,"Lower level",5X,"Energy (cm-1)" &
        &,1X,"F0VED0 (meV fm-2)",1X,"F0VED1 (meV fm-4)")'
    headerfp2mevlsj = '(1X,a,5X,a,5X,"Energy (cm-1)" &
        &,1X,"F0VED0 (meV fm-2)",1X,"F0VED1 (meV fm-4)")'

    headereghz = '(1X,"Upper level",5X,"Lower level",5X,"Energy (cm-1)",6X,"MS (GHz)",9X,"FS (GHz)",10X,"IS (GHz)")'
    headereghzlsj = '(1X,a,5X,a,5X,"Energy (cm-1)",6X,"MS (GHz)",9X,"FS (GHz)",10X,"IS (GHz)")'
    headeremhz = '(1X,"Upper level",5X,"Lower level",5X,"Energy (cm-1)",6X,"MS (MHz)",9X,"FS (MHz)",10X,"IS (MHz)")'
    headeremhzlsj = '(1X,a,5X,a,5X,"Energy (cm-1)",6X,"MS (MHz)",9X,"FS (MHz)",10X,"IS (MHz)")'
    headeremev = '(1X,"Upper level",5X,"Lower level",5X,"Energy (cm-1)",6X,"MS (meV)",9X,"FS (meV)",10X,"IS (meV)")'
    headeremevlsj = '(1X,a,5X,a,5X,"Energy (cm-1)",6X,"MS (meV)",9X,"FS (meV)",10X,"IS (meV)")'

    ! format1a = '(i4,a13,i3,a12,f12.2,1P,4d16.7)'
    ! format1b = '(i4,a13,i3,a12,f12.2,1P,4d16.7)'
    ! format1c = '(i4,a13,i3,a12,f12.2,1P,2d16.7)'
    ! format2 = '(i4,a13,i3,a12,f12.2,1P,3d16.7)'

    format1a = "(i4,a13,i3,a12,f12.2,1P,4d18.7)"
    format1alsj = "(a,a10,a,a10,f12.2,1P,4d18.7)"
    format1b = "(i4,a13,i3,a12,f12.2,1P,4d18.7)"
    format1blsj = "(a,a10,a,a10,f12.2,1P,4d18.7)"
    format1b2 = "(i4,a13,i3,a12,f12.2,1P,4d18.7,0P,f18.7)"
    format1b2lsj = "(a,a10,a,a10,f12.2,1P,4d18.7,0P,f18.7)"

    format1c = "(i4,a13,i3,a12,f12.2,1P,2d18.7)"
    format1clsj = "(a,a10,a,a10,f12.2,1P,2d18.7)"
    format1c2 = "(i4,a13,i3,a12,f12.2,1P,2d18.7,0P,f18.7)"
    format1c2lsj = "(a,a10,a,a10,f12.2,1P,2d18.7,0P,f18.7)"

    format2 = "(i4,a13,i3,a12,f12.2,1P,3d18.7)"
    format2lsj = "(a,a10,a,a10,f12.2,1P,3d18.7)"

    write(*,*) "Computes line frequency isotope shift parameters and/or energies using"
    write(*,*) "output files from ris4"
    write(*,*)
    write(*,*) "Input files: both states' .mcdhf.isotope-shift or .ci.isotope-shift reports, isodata"
    write(*,*) "Output file: <state1>.<state2>.(c)fi"

    ! Read the requested units and options.
    write(*,*)
    write(*,*) "Default settings (y/n)?: "
    read(*,*) def
    if(def/="y".and.def/="n") then
        write(*,*) 'Your input must be "y" or "n" (case sensitive). Try again!'
        stop
    end if

    call read_state_pair(state1, state2, roles=[character(11) :: "Upper state", "Lower state"])

    write(*,*) "Resulting isotope shifts from CI calculations (y/n)?: "
    read(*,*) ci
    if(ci/="y".and.ci/="n") then
        write(*,*) 'Your input must be "y" or "n" (case sensitive). Try again!'
        stop
    end if

    write(*,*) "Have electronic factors been calculated (y/n)?: "
    read(*,*) efac
    if(efac/="y".and.efac/="n") then
        write(*,*) 'Your input must be "y" or "n" (case sensitive). Try again!'
        stop
    end if

    if(def=="y") then
        typeflag = "both"
        unit  = "GHz"
        rel = "y"
        soph = "n"
    else if(def=="n") then
        write(*,*) "Compute IS parameters (para), IS energies (ener) or both (both)?: "
        read(*,*) typeflag
        if(typeflag/="para".and.typeflag/="ener".and.typeflag/="both") then
            write(*,*) 'Your input must be "para", "ener" or "both" (case sensitive). Try again!'
            stop
        end if

        write(*,*) "Units (GHz, MHz or meV)?: "
        read(*,*) unit
        if(unit/="GHz".and.unit/="MHz".and.unit/="meV") then
            write(*,*) 'Your input must be "GHz", "MHz", or "meV" (case sensitive). Try again!'
            stop
        end if

        write(*,*) "Use relativistically corrected mass shift parameters (y/n)?: "
        read(*,*) rel
        if(rel/="y".and.rel/="n") then
            write(*,*) 'Your input must be "y" or "n" (case sensitive). Try again!'
            stop
        end if

        if(typeflag=="ener".or.typeflag=="both") then
            write(*,*) "Use sophisticated model for radial moments (y/n)?: "
            read(*,*) soph
            if(soph/="y".and.soph/="n") then
                write(*,*) 'Your input must be "y" or "n" (case sensitive). Try again!'
                stop
            end if
        else if(typeflag=="para") then
            soph = "n"
        end if
    end if

    if ((typeflag=="ener").or.(typeflag=="both")) then
        if(soph=="y") then
            write(*,*) "Data for isotope 1"
            write(*,*) "Enter mass(amu),rms radius, t, omega, b20, b40: "
            read(*,*) m1,trms_1,t_1,w_1,b20_1,b40_1
            write(*,*) "Data for isotope 2"
            write(*,*) "Enter mass(amu),rms radius, t, omega, b20, b40: "
            read(*,*) m2,trms_2,t_2,w_2,b20_2,b40_2
        else if(soph=="n") then
            write(*,*) "Mass of isotope 1 (in amu)?: "
            read(*,*) m1
            write(*,*) "Mass of isotope 2 (in amu)?: "
            read(*,*) m2
            write(*,*) "d<r^2> (isotope 1 - isotope 2) (in fm^2)?: "
            read(*,*) dr2
            radius_difference = dr2
            if (m1 < m2) radius_difference = -dr2
        end if
        if (.not. (m1 > 0.0_real64 .and. m2 > 0.0_real64)) error stop "FICAL: isotope masses must be positive"
        inverse_mass_difference = -abs(m2-m1)/(m1*m2)
    end if

    ! Resolve the two state reports and the transition output path.

    file5 = lsj_label_input_path(state1, ci == "y")
    file6 = lsj_label_input_path(state2, ci == "y")
    file1 = isotope_shift_input_path(state1, ci == "y")
    file2 = isotope_shift_input_path(state2, ci == "y")
    if(ci=="y") then
        file3 = trim(state1)//"."//trim(state2)//".cfi"
    else
        file3 = trim(state1)//"."//trim(state2)//".fi"
    end if

    ! Select the output units and the supported field-shift terms.

    if(rel=="y") then
        relflag = 1
    else
        relflag = 0
    end if

    if(efac=="y") then
        efacflag = 1
    else
        efacflag = 0
    end if

    ! Load the electronic factors for each state.

    lsjpresent = 0
    nveclsj1 = 0
    nveclsj2 = 0
    inquire(file=file5, exist=exist1)
    inquire(file=file6, exist=exist2)
    if(exist1.and.exist2) then
        lsjpresent = 1
        call read_lsj_report(file5,nveclsj1,levellsj1,spinlsj1,parlsj1,csflabel1,lmax1)
        call read_lsj_report(file6,nveclsj2,levellsj2,spinlsj2,parlsj2,csflabel2,lmax2)
        lmax = max(lmax1,lmax2)
        uplevel = "Upper level                                                   "
        lowlevel = "Lower level                                                   "

    end if

    call read_isotope_shift_report(file1,unit,efacflag,relflag,first_state_level_count,level1,spinpar1,energy1,nms1,sms1,ms1,fs1)

    call read_isotope_shift_report(file2,unit,efacflag,relflag,second_state_level_count,level2,spinpar2,energy2,nms2,sms2,ms2,fs2)

    call read_isotope_data("isodata",file3,z)

    ! Write the nuclear input data to the report.
    if(typeflag/="para") then
        write(isotope_shift_output_unit,"(a)") " INPUT ISOTOPE DATA"
        if(soph=="n") then
            write(isotope_shift_output_unit,"(a22,f10.6,a)") " Mass of isotope 1: ",m1," amu"
            write(isotope_shift_output_unit,"(a22,f10.6,a)") " Mass of isotope 2: ",m2," amu"
            write(isotope_shift_output_unit,"(a22,f10.6,a)") " d<r^2>: ",dr2," fm^2"
        else
            write(isotope_shift_output_unit,*) "               Isotope 1    Isotope 2"
            write(isotope_shift_output_unit,"(a16,f10.6,f13.6)") " Mass [amu]: ",m1,m2
            write(isotope_shift_output_unit,"(a16,f10.6,f13.6)") " r_rms [fm]: ",trms_1,trms_2
            write(isotope_shift_output_unit,"(a16,f10.6,f13.6)") " t [fm]: ",t_1,t_2
            write(isotope_shift_output_unit,"(a16,f10.6,f13.6)") " omega: ",w_1,w_2
            write(isotope_shift_output_unit,"(a16,f10.6,f13.6)") " b_20: ",b20_1,b20_2
            write(isotope_shift_output_unit,"(a16,f10.6,f13.6)") " b_40: ",b40_1,b40_2
        end if
    end if
    write(isotope_shift_output_unit,*)

    if(efacflag==1) then
        fshift1(:first_state_level_count) = fs1(5,:first_state_level_count)
        fshift2(:second_state_level_count) = fs2(5,:second_state_level_count)
        fshift11(:first_state_level_count) = fs1(6,:first_state_level_count)
        fshift22(:second_state_level_count) = fs2(6,:second_state_level_count)
    else
        if(unit=="GHz") then
            fs1(1,:first_state_level_count) = fs1(1,:first_state_level_count)*z*4.9210908105e-3_real64
            fs2(1,:second_state_level_count) = fs2(1,:second_state_level_count)*z*4.9210908105e-3_real64
        else if(unit=="MHz") then
            fs1(1,:first_state_level_count) = fs1(1,:first_state_level_count)*z*4.9210908105_real64
            fs2(1,:second_state_level_count) = fs2(1,:second_state_level_count)*z*4.9210908105_real64
        else if(unit=="meV") then
            fs1(1,:first_state_level_count) = fs1(1,:first_state_level_count)*z*4.9210908105e-3_real64/mev2ghz
            fs2(1,:second_state_level_count) = fs2(1,:second_state_level_count)*z*4.9210908105e-3_real64/mev2ghz
        end if
        fshift1(:first_state_level_count) = fs1(1,:first_state_level_count)
        fshift2(:second_state_level_count) = fs2(1,:second_state_level_count)
        fshift11(:first_state_level_count) = 0.0_real64
        fshift22(:second_state_level_count) = 0.0_real64
    end if

    if(soph=="y") then
        call fermi_nuclear_distribution(trms_1,t_1,w_1,b20_1,b40_1,0, radmom1)
        call fermi_nuclear_distribution(trms_2,t_2,w_2,b20_2,b40_2,0, radmom2)


        fieldshift1(:) = 0.0_real64
        do i=1,first_state_level_count
            do j=1,4
                fieldshift1(i) = fieldshift1(i) + fs1(j,i)*(radmom1(j)-radmom2(j))
            end do
        end do
        fieldshift2(:) = 0.0_real64
        do i=1,second_state_level_count
            do j=1,4
                fieldshift2(i) = fieldshift2(i) + fs2(j,i)*(radmom1(j)-radmom2(j))
            end do
        end do

        ! Write the nuclear radial moments to the report.
        if(typeflag/="para") then
            write(isotope_shift_output_unit,"(a)") " NUCLEAR RADIAL MOMENTS"
            write(isotope_shift_output_unit,*) "               Isotope 1    Isotope 2     Isotope 1 - Isotope 2"
            write(isotope_shift_output_unit,"(a16,1P,D11.5,D13.5,D15.5)") " <r^2> [fm^2]: ",radmom1(1),radmom2(1),radmom1(1) &
                & -radmom2(1)
            write(isotope_shift_output_unit,"(a16,1P,D11.5,D13.5,D15.5)") " <r^4> [fm^4]: ",radmom1(2),radmom2(2),radmom1(2) &
                & -radmom2(2)
            write(isotope_shift_output_unit,"(a16,1P,D11.5,D13.5,D15.5)") " <r^6> [fm^6]: ",radmom1(3),radmom2(3),radmom1(3) &
                & -radmom2(3)
            write(isotope_shift_output_unit,"(a16,1P,D11.5,D13.5,D15.5)") " <r^8> [fm^8]: ",radmom1(4),radmom2(4),radmom1(4) &
                & -radmom2(4)
            write(isotope_shift_output_unit,*)
        end if

    end if

    ! Subtract the state factors to obtain line isotope-shift parameters.

    if((typeflag=="para").or.(typeflag=="both")) then
        write(isotope_shift_output_unit,*) "LINE MASS SHIFT PARAMETERS"
        write(isotope_shift_output_unit,*)
        if(unit=="GHz") then
            nmsconst = 29.9792458_real64/1822.888485_real64
            if(lsjpresent==0) then
                write(isotope_shift_output_unit,headermpghz)
            else
                write(isotope_shift_output_unit,headermpghzlsj) uplevel(1:lmax+8),lowlevel(1:lmax+8)
            end if
        else if(unit=="MHz") then
            nmsconst = 29979.2458_real64/1822.888485_real64
            if(lsjpresent==0) then
                write(isotope_shift_output_unit,headermpmhz)
            else
                write(isotope_shift_output_unit,headermpmhzlsj) uplevel(1:lmax+8),lowlevel(1:lmax+8)
            end if
        else if(unit=="meV") then
            nmsconst = 0.12398419739_real64/1822.888485_real64
            if(lsjpresent==0) then
                write(isotope_shift_output_unit,headermpmev)
            else
                write(isotope_shift_output_unit,headermpmevlsj) uplevel(1:lmax+8),lowlevel(1:lmax+8)
            end if
        end if
        write(isotope_shift_output_unit,*)
        do i=1,first_state_level_count
            do ii=1,nveclsj1
                if((level1(i)==levellsj1(ii)).and.(spinpar1(i)(6:9)==spinlsj1(ii)).and.(spinpar1(i)(11:11)==parlsj1(ii))) then
                    csflabel1scal = csflabel1(ii)
                end if
            end do
            do j=1,second_state_level_count
                do jj=1,nveclsj2
                    if((level2(j)==levellsj2(jj)).and.(spinpar2(j)(6:9)==spinlsj2(jj)).and.(spinpar2(j)(11:11)==parlsj2(jj))) then
                        csflabel2scal = csflabel2(jj)
                    end if
                end do
                if(energy1(i)>energy2(j)) then
                    tenergy = (energy1(i)-energy2(j))*hartree_to_wavenumber_cm
                    nmslaw = -1.0_real64*tenergy*nmsconst
                    if(state1==state2.and.j>i) then
                        if(lsjpresent==0) then
                            write(isotope_shift_output_unit,format1a) level1(i), adjustl(trim(spinpar1(i))),level2(j), &
                                adjustl(trim(spinpar2(j))), &
                                tenergy,nmslaw,nms1(i)-nms2(j),sms1(i)-sms2(j),ms1(i)-ms2(j)
                        else
                            write(isotope_shift_output_unit,format1alsj) csflabel1scal(1:lmax+3),adjustl(trim(spinpar1(i))), &
                                csflabel2scal(1:lmax+3), &
                                adjustl(trim(spinpar2(j))),tenergy,nmslaw,nms1(i)-nms2(j),sms1(i)-sms2(j),ms1(i)-ms2(j)
                        end if
                    else if(state1/=state2) then
                        if(lsjpresent==0) then
                            write(isotope_shift_output_unit,format1a) level1(i), adjustl(trim(spinpar1(i))),level2(j), &
                                adjustl(trim(spinpar2(j))), &
                                tenergy,nmslaw,nms1(i)-nms2(j),sms1(i)-sms2(j),ms1(i)-ms2(j)
                        else
                            write(isotope_shift_output_unit,format1alsj) csflabel1scal(1:lmax+3),adjustl(trim(spinpar1(i))), &
                                csflabel2scal(1:lmax+3), &
                                adjustl(trim(spinpar2(j))),tenergy,nmslaw,nms1(i)-nms2(j),sms1(i)-sms2(j),ms1(i)-ms2(j)
                        end if
                    end if
                else
                    tenergy = (energy2(j)-energy1(i))*hartree_to_wavenumber_cm
                    nmslaw = -1.0_real64*tenergy*nmsconst
                    if(state1==state2.and.j>i) then
                        if(lsjpresent==0) then
                            write(isotope_shift_output_unit,format1a) level2(j), adjustl(trim(spinpar2(j))),level1(i), &
                                adjustl(trim(spinpar1(i))), &
                                tenergy,nmslaw,nms2(j)-nms1(i),sms2(j)-sms1(i),ms2(j)-ms1(i)
                        else
                            write(isotope_shift_output_unit,format1alsj) csflabel2scal(1:lmax+3),adjustl(trim(spinpar2(j))), &
                                csflabel1scal(1:lmax+3), &
                                adjustl(trim(spinpar1(i))),tenergy,nmslaw,nms2(j)-nms1(i),sms2(j)-sms1(i),ms2(j)-ms1(i)
                        end if
                    else if(state1/=state2) then
                        if(lsjpresent==0) then
                            write(isotope_shift_output_unit,format1a) level2(j), adjustl(trim(spinpar2(j))),level1(i), &
                                adjustl(trim(spinpar1(i))), &
                                tenergy,nmslaw,nms2(j)-nms1(i),sms2(j)-sms1(i),ms2(j)-ms1(i)
                        else
                            write(isotope_shift_output_unit,format1alsj) csflabel2scal(1:lmax+3),adjustl(trim(spinpar2(j))), &
                                csflabel1scal(1:lmax+3), &
                                adjustl(trim(spinpar1(i))),tenergy,nmslaw,nms2(j)-nms1(i),sms2(j)-sms1(i),ms2(j)-ms1(i)
                        end if
                    end if
                end if
            end do
        end do
        write(isotope_shift_output_unit,*)
        write(isotope_shift_output_unit,*) "LINE FIELD SHIFT PARAMETERS"
        write(isotope_shift_output_unit,*)
        if(unit=="GHz") then
            if(lsjpresent==0) then
                write(isotope_shift_output_unit,headerfpghz)
            else
                write(isotope_shift_output_unit,headerfpghzlsj) uplevel(1:lmax+8),lowlevel(1:lmax+8)
            end if
        else if(unit=="MHz") then
            if(lsjpresent==0) then
                write(isotope_shift_output_unit,headerfpmhz)
            else
                write(isotope_shift_output_unit,headerfpmhzlsj) uplevel(1:lmax+8),lowlevel(1:lmax+8)
            end if
        else if(unit=="meV") then
            if(lsjpresent==0) then
                write(isotope_shift_output_unit,headerfpmev)
            else
                write(isotope_shift_output_unit,headerfpmevlsj) uplevel(1:lmax+8),lowlevel(1:lmax+8)
            end if
        end if
        write(isotope_shift_output_unit,*)
        do i=1,first_state_level_count
            do ii=1,nveclsj1
                if((level1(i)==levellsj1(ii)).and.(spinpar1(i)(6:9)==spinlsj1(ii)).and.(spinpar1(i)(11:11)==parlsj1(ii))) then
                    csflabel1scal = csflabel1(ii)
                end if
            end do
            do j=1,second_state_level_count
                do jj=1,nveclsj2
                    if((level2(j)==levellsj2(jj)).and.(spinpar2(j)(6:9)==spinlsj2(jj)).and.(spinpar2(j)(11:11)==parlsj2(jj))) then
                        csflabel2scal = csflabel2(jj)
                    end if
                end do
                if(energy1(i)>energy2(j)) then
                    tenergy = (energy1(i)-energy2(j))*hartree_to_wavenumber_cm
                    f0 = fs1(1,i)-fs2(1,j)
                    f2 = fs1(2,i)-fs2(2,j)
                    f4 = fs1(3,i)-fs2(3,j)
                    f6 = fs1(4,i)-fs2(4,j)
                    fieldshift = fshift1(i)-fshift2(j)
                    fieldshift22 = fshift11(i)-fshift22(j)
                    if(i==1.and.j==1) then
                        f0res = f0
                        f2res = f2
                    end if
                    linfac = linearity_discrepancy(f0,f2,f0res,f2res)

                    if(state1==state2.and.j>i) then
                        if(lsjpresent==0) then
                            write(isotope_shift_output_unit,format1b) level1(i), adjustl(trim(spinpar1(i))),level2(j), &
                                adjustl(trim(spinpar2(j))), &
                                tenergy,f0,f2,f4,f6
                        else
                            write(isotope_shift_output_unit,format1blsj) csflabel1scal(1:lmax+3),adjustl(trim(spinpar1(i))), &
                                csflabel2scal(1:lmax+3), &
                                adjustl(trim(spinpar2(j))),tenergy,f0,f2,f4,f6
                        end if
                    else if(state1/=state2) then
                        if(lsjpresent==0) then
                            write(isotope_shift_output_unit,format1b2) level1(i), adjustl(trim(spinpar1(i))),level2(j), &
                                adjustl(trim(spinpar2(j))), &
                                tenergy,f0,f2,f4,f6,linfac
                        else
                            write(isotope_shift_output_unit,format1b2lsj) csflabel1scal(1:lmax+3),adjustl(trim(spinpar1(i))), &
                                csflabel2scal(1:lmax+3), &
                                adjustl(trim(spinpar2(j))),tenergy,f0,f2,f4,f6,linfac
                        end if
                    end if
                else
                    tenergy = (energy2(j)-energy1(i))*hartree_to_wavenumber_cm
                    f0 = fs2(1,j)-fs1(1,i)
                    f2 = fs2(2,j)-fs1(2,i)
                    f4 = fs2(3,j)-fs1(3,i)
                    f6 = fs2(4,j)-fs1(4,i)
                    fieldshift = fshift2(j)-fshift1(i)
                    fieldshift22 = fshift22(j)-fshift11(i)
                    if(i==1.and.j==1) then
                        f0res = f0
                        f2res = f2
                    end if
                    linfac = linearity_discrepancy(f0,f2,f0res,f2res)

                    if(state1==state2.and.j>i) then
                        if(lsjpresent==0) then
                            write(isotope_shift_output_unit,format1b) level2(j), adjustl(trim(spinpar2(j))),level1(i), &
                                adjustl(trim(spinpar1(i))), &
                                tenergy,f0,f2,f4,f6
                        else
                            write(isotope_shift_output_unit,format1blsj) csflabel2scal(1:lmax+3),adjustl(trim(spinpar2(j))), &
                                csflabel1scal(1:lmax+3), &
                                adjustl(trim(spinpar1(i))),tenergy,f0,f2,f4,f6
                        end if
                    else if(state1/=state2) then
                        if(lsjpresent==0) then
                            write(isotope_shift_output_unit,format1b2) level2(j), adjustl(trim(spinpar2(j))),level1(i), &
                                adjustl(trim(spinpar1(i))), &
                                tenergy,f0,f2,f4,f6,linfac
                        else
                            write(isotope_shift_output_unit,format1b2lsj) csflabel2scal(1:lmax+3),adjustl(trim(spinpar2(j))), &
                                csflabel1scal(1:lmax+3), &
                                adjustl(trim(spinpar1(i))),tenergy,f0,f2,f4,f6,linfac
                        end if

                    end if
                end if
            end do
        end do

        write(isotope_shift_output_unit,*)
        if(unit=="GHz") then
            if(lsjpresent==0) then
                write(isotope_shift_output_unit,headerfp2ghz)
            else
                write(isotope_shift_output_unit,headerfp2ghzlsj) uplevel(1:lmax+8),lowlevel(1:lmax+8)
            end if
        else if(unit=="MHz") then
            if(lsjpresent==0) then
                write(isotope_shift_output_unit,headerfp2mhz)
            else
                write(isotope_shift_output_unit,headerfp2mhzlsj) uplevel(1:lmax+8),lowlevel(1:lmax+8)
            end if
        else if(unit=="meV") then
            if(lsjpresent==0) then
                write(isotope_shift_output_unit,headerfp2mev)
            else
                write(isotope_shift_output_unit,headerfp2mevlsj) uplevel(1:lmax+8),lowlevel(1:lmax+8)
            end if
        end if
        write(isotope_shift_output_unit,*)
        do i=1,first_state_level_count
            do ii=1,nveclsj1
                if((level1(i)==levellsj1(ii)).and.(spinpar1(i)(6:9)==spinlsj1(ii)).and.(spinpar1(i)(11:11)==parlsj1(ii))) then
                    csflabel1scal = csflabel1(ii)
                end if
            end do
            do j=1,second_state_level_count
                do jj=1,nveclsj2
                    if((level2(j)==levellsj2(jj)).and.(spinpar2(j)(6:9)==spinlsj2(jj)).and.(spinpar2(j)(11:11)==parlsj2(jj))) then
                        csflabel2scal = csflabel2(jj)
                    end if
                end do
                if(energy1(i)>energy2(j)) then
                    tenergy = (energy1(i)-energy2(j))*hartree_to_wavenumber_cm
                    f0 = fs1(1,i)-fs2(1,j)
                    f2 = fs1(2,i)-fs2(2,j)
                    f4 = fs1(3,i)-fs2(3,j)
                    f6 = fs1(4,i)-fs2(4,j)
                    fieldshift = fshift1(i)-fshift2(j)
                    fieldshift22 = fshift11(i)-fshift22(j)
                    if(i==1.and.j==1) then
                        fieldshiftres = fieldshift
                        fieldshift22res = fieldshift22
                    end if
                    linfac2 = linearity_discrepancy(fieldshift,fieldshift22,fieldshiftres,fieldshift22res)

                    if(state1==state2.and.j>i) then
                        if(lsjpresent==0) then
                            write(isotope_shift_output_unit,format1c) level1(i), adjustl(trim(spinpar1(i))),level2(j), &
                                adjustl(trim(spinpar2(j))), &
                                tenergy,fieldshift,fieldshift22
                        else
                            write(isotope_shift_output_unit,format1clsj) csflabel1scal(1:lmax+3),adjustl(trim(spinpar1(i))), &
                                csflabel2scal(1:lmax+3), &
                                adjustl(trim(spinpar2(j))),tenergy,fieldshift,fieldshift22
                        end if
                    else if(state1/=state2) then
                        if(lsjpresent==0) then
                            write(isotope_shift_output_unit,format1c2) level1(i), adjustl(trim(spinpar1(i))),level2(j), &
                                adjustl(trim(spinpar2(j))), &
                                tenergy,fieldshift,fieldshift22,linfac2
                        else
                            write(isotope_shift_output_unit,format1c2lsj) csflabel1scal(1:lmax+3),adjustl(trim(spinpar1(i))), &
                                csflabel2scal(1:lmax+3), &
                                adjustl(trim(spinpar2(j))),tenergy,fieldshift,fieldshift22,linfac2
                        end if
                    end if
                else
                    tenergy = (energy2(j)-energy1(i))*hartree_to_wavenumber_cm
                    f0 = fs2(1,j)-fs1(1,i)
                    f2 = fs2(2,j)-fs1(2,i)
                    f4 = fs2(3,j)-fs1(3,i)
                    f6 = fs2(4,j)-fs1(4,i)
                    fieldshift = fshift2(j)-fshift1(i)
                    fieldshift22 = fshift22(j)-fshift11(i)
                    if(i==1.and.j==1) then
                        fieldshiftres = fieldshift
                        fieldshift22res = fieldshift22
                    end if
                    linfac2 = linearity_discrepancy(fieldshift,fieldshift22,fieldshiftres,fieldshift22res)

                    if(state1==state2.and.j>i) then
                        if(lsjpresent==0) then
                            write(isotope_shift_output_unit,format1c) level2(j), adjustl(trim(spinpar2(j))),level1(i), &
                                adjustl(trim(spinpar1(i))), &
                                tenergy,fieldshift,fieldshift22
                        else
                            write(isotope_shift_output_unit,format1clsj) csflabel2scal(1:lmax+3),adjustl(trim(spinpar2(j))), &
                                csflabel1scal(1:lmax+3), &
                                adjustl(trim(spinpar1(i))),tenergy,fieldshift,fieldshift22
                        end if
                    else if(state1/=state2) then
                        if(lsjpresent==0) then
                            write(isotope_shift_output_unit,format1c2) level2(j), adjustl(trim(spinpar2(j))),level1(i), &
                                adjustl(trim(spinpar1(i))), &
                                tenergy,fieldshift,fieldshift22,linfac2
                        else
                            write(isotope_shift_output_unit,format1c2lsj) csflabel2scal(1:lmax+3),adjustl(trim(spinpar2(j))), &
                                csflabel1scal(1:lmax+3), &
                                adjustl(trim(spinpar1(i))),tenergy,fieldshift,fieldshift22,linfac2
                        end if
                    end if
                end if
            end do
        end do

    end if

    ! Combine the line factors with the nuclear data to obtain isotope shifts.

    write(isotope_shift_output_unit,*)
    if((typeflag=="ener").or.(typeflag=="both")) then
        write(isotope_shift_output_unit,*) "LINE ISOTOPE SHIFT ENERGIES"
        write(isotope_shift_output_unit,*)
        if(unit=="GHz") then
            if(lsjpresent==0) then
                write(isotope_shift_output_unit,headereghz)
            else
                write(isotope_shift_output_unit,headereghzlsj) uplevel(1:lmax+8),lowlevel(1:lmax+8)
            end if
        else if(unit=="MHz") then
            if(lsjpresent==0) then
                write(isotope_shift_output_unit,headeremhz)
            else
                write(isotope_shift_output_unit,headeremhzlsj) uplevel(1:lmax+8),lowlevel(1:lmax+8)
            end if
        else if(unit=="meV") then
            if(lsjpresent==0) then
                write(isotope_shift_output_unit,headeremev)
            else
                write(isotope_shift_output_unit,headeremevlsj) uplevel(1:lmax+8),lowlevel(1:lmax+8)
            end if
        end if

        write(isotope_shift_output_unit,*)
        do i=1,first_state_level_count
            do ii=1,nveclsj1
                if((level1(i)==levellsj1(ii)).and.(spinpar1(i)(6:9)==spinlsj1(ii)).and.(spinpar1(i)(11:11)==parlsj1(ii))) then
                    csflabel1scal = csflabel1(ii)
                end if
            end do
            do j=1,second_state_level_count
                do jj=1,nveclsj2
                    if((level2(j)==levellsj2(jj)).and.(spinpar2(j)(6:9)==spinlsj2(jj)).and.(spinpar2(j)(11:11)==parlsj2(jj))) then
                        csflabel2scal = csflabel2(jj)
                    end if
                end do
                ! The line uses upper minus lower; isotope shifts use heavy minus light.
                masshift = (ms1(i)-ms2(j))*inverse_mass_difference
                if (soph == "y") then
                    fieldshift = fieldshift1(i)-fieldshift2(j)
                    if (m1 < m2) fieldshift = -fieldshift
                else
                    fieldshift = (fshift1(i)-fshift2(j))*radius_difference + &
                        (fshift11(i)-fshift22(j))*radius_difference**2
                end if
                if (energy1(i) <= energy2(j)) then
                    masshift = -masshift
                    fieldshift = -fieldshift
                end if
                isotopeshift = masshift + fieldshift
                if(energy1(i)>energy2(j)) then
                    tenergy = (energy1(i)-energy2(j))*hartree_to_wavenumber_cm
                    if(state1==state2.and.j>i) then
                        if(lsjpresent==0) then
                            write(isotope_shift_output_unit,format2) level1(i),adjustl(trim(spinpar1(i))),level2(j), &
                                adjustl(trim(spinpar2(j))),tenergy, &
                                masshift,fieldshift,isotopeshift
                        else
                            write(isotope_shift_output_unit,format2lsj) csflabel1scal(1:lmax+3),adjustl(trim(spinpar1(i))), &
                                csflabel2scal(1:lmax+3), &
                                adjustl(trim(spinpar2(j))),tenergy,masshift,fieldshift,isotopeshift
                        end if
                    else if(state1/=state2) then
                        if(lsjpresent==0) then
                            write(isotope_shift_output_unit,format2) level1(i),adjustl(trim(spinpar1(i))),level2(j), &
                                adjustl(trim(spinpar2(j))),tenergy, &
                                masshift,fieldshift,isotopeshift
                        else
                            write(isotope_shift_output_unit,format2lsj) csflabel1scal(1:lmax+3),adjustl(trim(spinpar1(i))), &
                                csflabel2scal(1:lmax+3), &
                                adjustl(trim(spinpar2(j))),tenergy,masshift,fieldshift,isotopeshift
                        end if
                    end if
                else
                    tenergy = (energy2(j)-energy1(i))*hartree_to_wavenumber_cm
                    if(state1==state2.and.j>i) then
                        if(lsjpresent==0) then
                            write(isotope_shift_output_unit,format2) level2(j),adjustl(trim(spinpar2(j))),level1(i), &
                                adjustl(trim(spinpar1(i))),tenergy, &
                                masshift,fieldshift,isotopeshift
                        else
                            write(isotope_shift_output_unit,format2lsj) csflabel2scal(1:lmax+3),adjustl(trim(spinpar2(j))), &
                                csflabel1scal(1:lmax+3), &
                                adjustl(trim(spinpar1(i))),tenergy,masshift,fieldshift,isotopeshift
                        end if
                    else if(state1/=state2) then
                        if(lsjpresent==0) then
                            write(isotope_shift_output_unit,format2) level2(j),adjustl(trim(spinpar2(j))),level1(i), &
                                adjustl(trim(spinpar1(i))),tenergy, &
                                masshift,fieldshift,isotopeshift
                        else
                            write(isotope_shift_output_unit,format2lsj) csflabel2scal(1:lmax+3),adjustl(trim(spinpar2(j))), &
                                csflabel1scal(1:lmax+3), &
                                adjustl(trim(spinpar1(i))),tenergy,masshift,fieldshift,isotopeshift
                        end if
                    end if
                end if
            end do
        end do
    end if
    write(isotope_shift_output_unit,*)
    close(isotope_shift_output_unit)
    write(*,*)
    write(*,*) "program FICAL finished ..."
    write(*,*)
    write(*,*) "Isotope shift parameters/energies written to file ", file3
    call finish_run_timer (start_clock_count, "FICAL")
contains

    pure function linearity_discrepancy(first,second,ref_first,ref_second) result(value)
        ! Module dependencies
        use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
        implicit none (type, external)

        ! Arguments
        real(real64), intent(in) :: first, second, ref_first, ref_second
        real(real64) :: value, denominator

        denominator = second*ref_first
        if (abs(denominator) > 0.0_real64) then
            value = 1000.0_real64*(first*ref_second/denominator-1.0_real64)
        else
            ! Electron-density-only reports have zero higher-order coefficients.
            value = ieee_value(0.0_real64,ieee_quiet_nan)
        end if
    end function linearity_discrepancy

    subroutine read_isotope_data(name,nameout,z)
        implicit none (type, external)

        ! Arguments and local declarations
        integer :: isotope_shift_input_unit
        character(*), intent(in) :: name, nameout
        real(real64), intent(out) :: z
        real(real64) :: nuclear_mass_number, a,c,radmom3(4)

        open (newunit=isotope_shift_input_unit,file=name,status="old",action="read",form="formatted")
        call read_isotope_value(isotope_shift_input_unit,name,"atomic number",z)
        call read_isotope_value(isotope_shift_input_unit,name,"mass number",nuclear_mass_number)
        call read_isotope_value(isotope_shift_input_unit,name,"Fermi diffuseness",a)
        call read_isotope_value(isotope_shift_input_unit,name,"Fermi radius",c)
        if (nuclear_mass_number == 0.0_real64) then
            a = 0.0_real64
            c = 0.0_real64
            radmom3 = 0.0_real64
        else
            if (a <= 0.0_real64 .or. c <= 0.0_real64) error stop "FICAL: invalid Fermi parameters in "//trim(name)
            a = a*4.0_real64*log(3.0_real64)
            call fermi_nuclear_distribution(c,a,0.0_real64,0.0_real64,0.0_real64,1,radmom3)
        end if
        open (newunit=isotope_shift_output_unit,file=nameout,status="replace",form="formatted")
        write(isotope_shift_output_unit,*)
        write(isotope_shift_output_unit,"(a)") " REFERENCE ISOTOPE DATA FROM ISODATA"
        write(isotope_shift_output_unit,"(a,f10.6)") " Atomic number: ", z
        if (nuclear_mass_number /= 0.0_real64) then
            write(isotope_shift_output_unit,"(a)") " Fermi nucleus: "
        else
            write(isotope_shift_output_unit,"(a)") " Point nucleus: "
        end if
        write(isotope_shift_output_unit,"(a16,f10.6,a)") "c: ", c, " fm"
        write(isotope_shift_output_unit,"(a16,f10.6,a)") "r_rms: ", sqrt(radmom3(1)), " fm"
        write(isotope_shift_output_unit,"(a16,f10.6,a)") "t: ", a, " fm"
        write(isotope_shift_output_unit,*)

        close(isotope_shift_input_unit)

    end subroutine read_isotope_data

    subroutine read_isotope_value(input,name,section,value)
        integer, intent(in) :: input
        character(*), intent(in) :: name, section
        real(real64), intent(out) :: value
        character(200) :: record
        integer :: status

        call next_report_line(input,name,section,record)
        call next_report_line(input,name,section,record)
        read(record,*,iostat=status) value
        if (status /= 0) error stop "FICAL: invalid "//section//" in "//trim(name)
        if (.not.ieee_is_finite(value)) error stop "FICAL: nonfinite "//section//" in "//trim(name)
    end subroutine read_isotope_value

    subroutine next_report_line(input,name,section,record)
        integer, intent(in) :: input
        character(*), intent(in) :: name, section
        character(*), intent(out) :: record
        integer :: status

        do
            read(input,"(a)",iostat=status) record
            if (status /= 0) error stop "FICAL: reading "//section//" in "//trim(name)
            if (len_trim(record) /= 0) return
        end do
    end subroutine next_report_line

    subroutine expect_heading(input,name,heading)
        integer, intent(in) :: input
        character(*), intent(in) :: name, heading
        character(200) :: record

        call next_report_line(input,name,heading,record)
        if (index(record,heading) == 0) error stop "FICAL: expected "//heading//" in "//trim(name)
    end subroutine expect_heading

    subroutine expect_field_units(input,name,labels)
        integer, intent(in) :: input
        character(*), intent(in) :: name, labels(:)
        character(200) :: record
        integer :: column, position

        call next_report_line(input,name,"field factor units",record)
        record = adjustl(record)
        ! Older abbreviated headings specify the first column only.
        if (record == labels(1)) return
        do column=1,size(labels)
            position = index(record,trim(labels(column)))
            if (position /= 1) error stop "FICAL: expected "//trim(labels(column))//" in "//trim(name)
            record = adjustl(record(len_trim(labels(column))+1:))
        end do
    end subroutine expect_field_units

    subroutine read_state_row(input,name,section,values,root,symmetry,atomic_units)
        integer, intent(in) :: input
        character(*), intent(in) :: name, section
        real(real64), intent(out) :: values(:)
        integer, intent(out) :: root
        character(11), intent(out) :: symmetry
        logical, intent(in) :: atomic_units
        character(200) :: record
        integer :: field, start, status

        call next_report_line(input,name,section,record)
        read(record,"(i4,a11)",iostat=status) root,symmetry
        if (status /= 0 .or. len_trim(record(:4)) == 0) error stop "FICAL: invalid level in "//section//" in "//trim(name)
        if (len_trim(record) < 20+20*size(values)) error stop "FICAL: incomplete "//section//" in "//trim(name)
        do field=1,size(values)
            start = 21+20*(field-1)
            if (record(start:start+19) == "") error stop "FICAL: missing "//section//" value in "//trim(name)
            read(record(start:start+19),"(d20.10)",iostat=status) values(field)
            if (status /= 0) error stop "FICAL: invalid "//section//" in "//trim(name)
            if (.not.ieee_is_finite(values(field))) error stop "FICAL: nonfinite "//section//" in "//trim(name)
        end do
        if (atomic_units .and. adjustl(record(21+20*size(values):)) /= "(a.u.)") &
            error stop "FICAL: expected atomic units for "//section//" in "//trim(name)
    end subroutine read_state_row

    subroutine read_isotope_shift_report(name,unit,eflag,rflag,nvec,level,spinpar,energy,nms,sms,ms,fs)
        implicit none (type, external)

        ! Arguments
        character(*), intent(in) :: name, unit
        integer, intent(in) :: eflag, rflag
        integer, intent(out) :: nvec, level(:)
        character(len=11), intent(out) :: spinpar(:)
        real(real64), intent(out) :: energy(:), nms(:), sms(:), ms(:), fs(:,:)

        ! Local declarations
        integer :: input, i, component, operator, root, status
        character(11) :: symmetry
        character(3), parameter :: mass_names(2) = ["NMS","SMS"]

        real(real64) :: ghz_to_output, terms(3)
        character(len=200) :: string

        select case (unit)
        case ("GHz")
            ghz_to_output = 1.0_real64
        case ("MHz")
            ghz_to_output = 1000.0_real64
        case ("meV")
            ghz_to_output = 1.0_real64/mev2ghz
        case default
            error stop "FICAL: unsupported output unit "//trim(unit)
        end select
        component = 1
        if (rflag == 1) component = 3

        open(newunit=input, file=name, status="old", action="read", form="formatted")
        call next_report_line(input,name,"eigenvalue count",string)
        if (index(string,"Number of eigenvalues:") == 0) error stop "FICAL: missing eigenvalue count in "//trim(name)
        read(string,"(24x,i3)",iostat=status) nvec
        if (status /= 0) error stop "FICAL: invalid eigenvalue count in "//trim(name)
        if (nvec < 1 .or. nvec > size(level)) error stop "FICAL: invalid eigenvalue count in "//trim(name)
        call expect_heading(input,name,"Energy")
        do i=1,nvec
            call read_state_row(input,name,"energy",energy(i:i),level(i),spinpar(i),.true.)
        end do
        do operator=1,2
            if (operator == 1) then
                call expect_heading(input,name,"Normal mass shift parameter")
            else
                call expect_heading(input,name,"Specific mass shift parameter")
            end if
            do i=1,nvec
                call expect_heading(input,name,"<K^1>")
                call read_state_row(input,name,mass_names(operator),terms,root,symmetry,.true.)
                if (root /= level(i) .or. symmetry /= spinpar(i)) error stop "FICAL: mismatched "//mass_names(operator)// &
                    " level in "//trim(name)
                call read_mass_shift(input,name,mass_names(operator),terms)
                if (operator == 1) then
                    nms(i) = terms(component)*ghz_to_output
                else
                    sms(i) = terms(component)*ghz_to_output
                end if
            end do
        end do
        ms(:nvec) = nms(:nvec)+sms(:nvec)
        call expect_heading(input,name,"Electron density in atomic units")
        call expect_heading(input,name,"(a.u.)")
        do i=1,nvec
            call read_state_row(input,name,"contact density",fs(1:1,i),root,symmetry,.false.)
            if (root /= level(i) .or. symmetry /= spinpar(i)) error stop "FICAL: mismatched density level in "//trim(name)
            fs(2:6,i) = 0.0_real64
        end do
        if (eflag == 1) then
            call expect_heading(input,name,"electronic factors")
            call expect_field_units(input,name,["F0 (GHz/fm^2)","F2 (GHz/fm^4)","F4 (GHz/fm^6)","F6 (GHz/fm^8)"])
            do i=1,nvec
                call read_state_row(input,name,"field factors",fs(1:4,i),root,symmetry,.false.)
                if (root /= level(i) .or. symmetry /= spinpar(i)) error stop "FICAL: mismatched field level in "//trim(name)
            end do
            call expect_heading(input,name,"electronic factors")
            call expect_field_units(input,name,["F0VED0 (GHz/fm^2)","F0VED1 (GHz/fm^4)"])
            do i=1,nvec
                call read_state_row(input,name,"corrected field factors",fs(5:6,i),root,symmetry,.false.)
                if (root /= level(i) .or. symmetry /= spinpar(i)) error stop "FICAL: mismatched corrected level in "//trim(name)
            end do
            fs(:,:nvec) = fs(:,:nvec)*ghz_to_output
        end if

        close(input)
    end subroutine read_isotope_shift_report

    subroutine read_mass_shift(input,name,parameter_name,terms)
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: input
        character(*), intent(in) :: name, parameter_name
        real(real64), intent(out) :: terms(3)

        ! Local declarations
        character(len=200) :: record
        character(len=200) :: message
        integer :: status

        ! Formats 315 and 316 have the same three D20.10 fields, followed by
        ! a unit label. Read the label on every row, since even levels in one
        ! RIS4 file can use different units.
        read(input,"(a)",iostat=status,iomsg=message) record
        if (status /= 0) error stop "FICAL: reading "//parameter_name//" in "//trim(name)//": "//trim(message)
        if (record(21:40) == "" .or. record(41:60) == "" .or. record(61:80) == "") &
            error stop "FICAL: incomplete "//parameter_name//" record in "//trim(name)
        read(record,"(20x,3d20.10)",iostat=status,iomsg=message) terms
        if (status /= 0) error stop "FICAL: invalid "//parameter_name//" record in "//trim(name)//": "//trim(message)

        if (.not.all(ieee_is_finite(terms))) error stop "FICAL: nonfinite "//parameter_name//" in "//trim(name)
        select case (adjustl(record(81:)))
        case ("(MHz u)")
            terms = terms/1000.0_real64
        case ("(GHz u)")
        ! Already in GHz u.
        case default
            error stop "FICAL: unsupported "//parameter_name//" unit in "//trim(name)//" ["//trim(adjustl(record(81:)))//"]"
        end select
    end subroutine read_mass_shift

    subroutine fermi_nuclear_distribution(trms,t,w,b20,b40,flag, radmom)
        implicit none (type, external)

        ! Local declarations
        integer :: m, m2, i, j, k, flag
        real(real64) :: c0, a, t, w, b20, b40
        real(real64) :: rmax, h, h2, rms, trms, drms

        ! Local constants
        real(real64), parameter :: pi = 3.14159265359_real64

        ! Local declarations
        real(real64) :: dnorm, norm, dinte, inte
        real(real64) :: dr2, r2, dr4, r4, dr6, r6, dr8, r8, radmom(4)
        real(real64) :: cv(0:2000), th(0:2000)
        real(real64), allocatable :: r(:), rho(:)

        a = t/(4.0_real64*log(3.0_real64))

        rmax = 100.0_real64
        ! h = 0.1d0
        h = 0.01_real64
        m = int(rmax/h)
        m = m + mod(m,2)

        allocate ( r(0:m), rho(0:m) )

        c0 = trms
        k = 1
        drms = 1.0_real64
        do while (drms >= 1.0e-8_real64)
            if(k>1) then
                c0 = trms/rms*c0
            end if
            do i=0, m
                r(i) = h*i
                if(abs(b20)<0.001_real64.and.abs(b40)<0.001_real64) then
                    rho(i) = 4.0_real64*pi*(1.0_real64 + w*r(i)**2.0_real64/c0**2.0_real64)/(1.0_real64+exp((r(i)-c0)/a))
                end if
            end do

            if(abs(b20) > 0.0_real64.or.abs(b40) > 0.0_real64) then
                m2 = 360
                h2 = pi/360.0_real64
                do i=0, m2, 1
                    th(i) = i*pi/360.0_real64
                    cv(i) = c0*(1.0_real64 + b20*sqrt(5.0_real64/(16.0_real64*pi))*(3.0_real64*(cos(th(i)))**2.0_real64 - &
                        & 1.0_real64) &
                        + b40*3.0_real64/(16.0_real64*sqrt(pi))*(35.0_real64*(cos(th(i)))**4.0_real64 - &
                        30.0_real64*(cos(th(i)))**2.0_real64 + 3.0_real64))
                end do

                do i=0, m, 1
                    inte = 0.0
                    do j=1, m2-1, 2
                        dinte = h2/3.0_real64 &
                            *((1.0_real64 + w*r(i)**2.0_real64/c0**2.0_real64)/(1.0_real64 + exp((r(i) - cv(j-1))/a))*sin(th(j-1)) &
                            + 4.0_real64*(1.0_real64 + w*r(i)**2.0_real64/c0**2.0_real64)/(1.0_real64 + exp((r(i) - cv(j))/a)) &
                            & *sin(th(j)) &
                            + (1.0_real64 + w*r(i)**2.0_real64/c0**2.0_real64)/(1.0_real64 + exp((r(i) - cv(j+1))/a))*sin(th(j+1)))
                        inte = inte + dinte
                    end do
                    rho(i) = 2.0_real64*pi*inte
                end do
            end if

            norm = 0.0_real64
            r2 = 0.0_real64
            r4 = 0.0_real64
            r6 = 0.0_real64
            r8 = 0.0_real64
            do i=1, m-1, 2
                dnorm = (h/3.0_real64)*(rho(i-1)*r(i-1)**2.0_real64 + 4.0_real64*rho(i)*r(i)**2.0_real64 + rho(i+1)*r(i+1) &
                    & **2.0_real64)
                dr2 = (h/3.0_real64)*(rho(i-1)*r(i-1)**4.0_real64 + 4.0_real64*rho(i)*r(i)**4.0_real64 + rho(i+1)*r(i+1) &
                    & **4.0_real64)
                dr4 = (h/3.0_real64)*(rho(i-1)*r(i-1)**6.0_real64 + 4.0_real64*rho(i)*r(i)**6.0_real64 + rho(i+1)*r(i+1) &
                    & **6.0_real64)
                dr6 = (h/3.0_real64)*(rho(i-1)*r(i-1)**8.0_real64 + 4.0_real64*rho(i)*r(i)**8.0_real64 + rho(i+1)*r(i+1) &
                    & **8.0_real64)
                dr8 = (h/3.0_real64)*(rho(i-1)*r(i-1)**10.0_real64 + 4.0_real64*rho(i)*r(i)**10.0_real64 + rho(i+1)*r(i+1) &
                    & **10.0_real64)
                norm = norm + dnorm
                r2 = r2 + dr2
                r4 = r4 + dr4
                r6 = r6 + dr6
                r8 = r8 + dr8
            end do
            radmom(1) = r2/norm
            rms = sqrt(radmom(1))
            radmom(2) = r4/norm
            radmom(3) = r6/norm
            radmom(4) = r8/norm
            norm = 1.0_real64/norm
            k = k + 1
            if(flag==0) then
                drms = abs(trms-rms)
            else
                drms = 0.0_real64
            end if
        end do

        deallocate ( r, rho )

        return

    end subroutine fermi_nuclear_distribution

    subroutine read_lsj_report(name,nvec,level,spin,par,csflabel,lmax)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: iostat_end
        implicit none (type, external)

        ! Arguments and local declarations
        integer :: lsj_label_input_unit
        integer :: i, nvec,level(500),lmax,io_status
        character(len=256) :: io_message
        character(len=4) :: spin(500)
        character(len=1) :: par(500)
        character(len=10) :: spinpar(500)
        character(len=36) :: dummy

        character(*), intent(in) :: name
        character(len=200) :: string
        character(len=200) :: csflabel(500)

        open (newunit=lsj_label_input_unit,file=name,status="unknown",form="formatted")

        lmax = 0
        read(lsj_label_input_unit,"(a)") string
        i=1
        do
            read(lsj_label_input_unit,"(i4,a10)",iostat=io_status,iomsg=io_message) level(i),spinpar(i)
            if (io_status == iostat_end) exit
            if (io_status /= 0) error stop trim(io_message)
            if(spinpar(i)(10:10)=="+".or.spinpar(i)(10:10)=="-") then
                spin(i) = spinpar(i)(1:4)
                par(i) = spinpar(i)(10:10)
                read(lsj_label_input_unit,"(a36,a)") dummy,csflabel(i)
                if(len(trim(csflabel(i)))>lmax) lmax = len(trim(csflabel(i)))
                i = i + 1
                !    spin(i) = spinpar(i)(1:4)
                !    par(i) = spinpar(i)(10:10)
            end if

            ! i = i + 1
        end do
        close(lsj_label_input_unit)
        nvec = i-1
    end subroutine read_lsj_report

end program isotope_shift_transition_main
