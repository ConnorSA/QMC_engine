!> @brief Functions and subroutines for implementing Markov chain Monte Carlo
module mcmc
    use shared_constants
    implicit none
    !MCMC only constants go here, try using private to keep namespace clean

    !going to try using mainly local(in function) variables to (hopefully) make
    !parallelisation cleaner/easier, despite what my OO brain wants to do
    !Alisdair

    !required inputs:log_rho(log of square magnitude of waveform), x_0, n_steps, n_burned, thinning_interval, s(step size), seed
    !output: array of samples from rho
    !to find integral evaluate local energy at all points in output, sum, normlise by length of output

            !decribes shape of input functions
    abstract interface
    function log_rho_interface(x, dof) result(out)
        !log target dist
        use shared_constants
        real(dp), dimension(:), intent(in) :: x
        real(dp), dimension(:), intent(in) :: dof
        real(dp) :: out
    end function log_rho_interface
    end interface

    private :: mcmc_propose, mcmc_log_prop_dist, mcmc_accept

    contains


    !> @brief A Gaussian random walk proposal.
    !! Uses Box-Mueller to generate random draws
    !! assumes PRNG is already intialised
    !ensure PRNG is initialised before this is used
    function mcmc_propose(x,s,density_dimension) result(x_prop)
        implicit none
        integer, intent(in) :: density_dimension !> dimension of measure to be sampled from, same dimension of x
        real(dp), dimension(density_dimension), intent(in) :: x !> current point/mean of normal
        real(dp), intent(in) :: s  !> step size parameter/varience of normal
        real(dp), dimension(2*((density_dimension+1)/2)) :: u !uniform random sample,
        !note the dimension is the nearest even number above the given dimension
        !integer division will autmatically round towards 0
        real(dp), dimension(density_dimension) :: xi !normal random sample
        real(dp), dimension(density_dimension) :: x_prop !> proposal/draw from N(x,s)
        integer :: i
        !proposes, uses same dist as described in log_prop_dist
        !parametrising ititially with just a step size,
        !if proposal is changed, might need updated

        !implementing box muller to generate a normal
        call random_number(u)
        do i=1,density_dimension/2
            xi(i) = sqrt(-2.0_dp*log(u(i)))*cos(2.0_dp*pi*u(i*2))
            xi(i*2) = sqrt(-2.0_dp*log(u(i)))*sin(2.0_dp*pi*u(i*2))
        end do
        if (modulo(density_dimension,2).ne.0) then
            xi(density_dimension) = sqrt(-2.0_dp*log(u(density_dimension)))*cos(2.0_dp*pi*u(density_dimension+1))
        end if

        !basic GRW
        x_prop = x + s*xi

    end function mcmc_propose

    !> @brief log likelihood of sampling y given x, with mcmc_propose
    pure function mcmc_log_prop_dist(x,y,s) result(out)
        !log liklihood x proposed at y
        implicit none
        real(dp), dimension(:), intent(in) :: x
        real(dp), dimension(:), intent(in) :: y
        real(dp), intent(in) :: s
        real(dp) :: out

        !unnormalised, as not needed, doing 1D step not nessesearily best but is easy
        out= -norm2(x-y)**2*(1/s)

    end function mcmc_log_prop_dist

    !> @brief Acceptance probaility of x_prop, given x, using Metropolis-Hastings ratio
    subroutine mcmc_accept(a,x,x_prop,log_rho,s,rho_x,rho_x_prop,dof_coefficients)
        !computes accpect prob in [0,1]
        implicit none
        procedure(log_rho_interface) :: log_rho !>the (log) density to be sampled from, takes input, x, and parametrisation dof_coefficients
        real(dp), dimension(:), intent(in) :: x, x_prop !>the current and proposed states
        real(dp), optional, intent(in) :: rho_x !>log_rho(x) can be saved between runs to reduce calls
        real(dp), optional, intent(out) :: rho_x_prop !>log_rho(x_prop)
        real(dp), dimension(:), intent(in) :: dof_coefficients !> the parameterisation of log_rho
        real(dp), intent(in) :: s !>stepsize parameter, contols log_prop_dist
        real(dp), intent(out) :: a !>$\mathbb{P}[\textrm{accept} x_{prop}|x]$
        real(dp) :: rho_x_, rho_x_prop_

        if (present(rho_x)) then
            rho_x_ = rho_x
        else
            rho_x_ = log_rho(x, dof_coefficients)
        end if

        rho_x_prop_ = log_rho(x_prop, dof_coefficients)

        if (present(rho_x_prop)) rho_x_prop=rho_x_prop_

        a = (mcmc_log_prop_dist(x,x_prop,s)+rho_x_prop_)-(mcmc_log_prop_dist(x_prop,x,s)+rho_x_)

    end subroutine mcmc_accept

    !docs say this should initialise a new PRNG for each thread,
    !should be ok for parallelisation
    !> @brief The main routine for MCMC, generates (n_steps-n_burned)/thinning_interval samples from a distribution rho
    subroutine mcmc_sample(samples, log_rho, x_0, n_steps, n_burned, thinning_interval,s, e_code,&
         dof_coefficients, density_dimension, average_accept, seed)
        implicit none
        procedure(log_rho_interface) :: log_rho !>the (log) density to be sampled from, takes input, x, and parametrisation dof_coefficients
        integer, intent(in) :: density_dimension !>size of x
        integer, intent(in) :: n_steps !>total number of steps to take
        integer, intent(in) :: n_burned !>number of initial steps to ignore, slightly reduces convergence rate but saves a lot of memory
        integer, intent(in) :: thinning_interval !>interval between saving steps, decreases correlation between saved samples, and saves memory
        integer, intent(in), dimension(:), optional :: seed !>the random seed, will use current state if not set
        real(dp), dimension(:), intent(in) :: dof_coefficients!> the parameterisation of log_rho
        real(dp), intent(in) :: s !>stepsize parameter, contols log_prop_dist and mcmc_propose
        real(dp), dimension(density_dimension), intent(in) :: x_0 !>a starting point for the chain
        integer, intent(out) :: e_code !>an error code, 0 is successs, 1, is a NaN accpetance, 2 is a negative acceptance, 4 is array size mismatch
        real(dp), intent(out), optional :: average_accept !>average rate of accpetance of proposals, optimal is ~0.24
        real(dp), dimension((n_steps-n_burned)/thinning_interval+1,density_dimension), intent(out) :: samples !>output, an array conatining all sampled points
        real(dp), dimension(density_dimension) :: x,x_prop
        real(dp) :: a,a_ave,rand_sample,log_rho_x,log_rho_x_prop
        integer :: loop_index,out_array_index,thinning_counter

        !e_code 0 => no error
        e_code=0

        x = x_0
        out_array_index = 0
        thinning_counter = 0
        samples(1,:) = x_0
        a_ave = 0.0_dp
        log_rho_x = log_rho(x_0,dof_coefficients)

        !uses current random state if none is given
        if (present(seed)) then
            call random_seed(put=seed)
        end if

        do loop_index=1,n_steps
            x_prop = mcmc_propose(x,s,density_dimension)
            call mcmc_accept(a, x, x_prop, log_rho, s,log_rho_x,log_rho_x_prop, dof_coefficients)
            a = exp(a)
            call random_number(rand_sample)
            if (a>rand_sample) then
                x = x_prop
                log_rho_x = log_rho_x_prop !using datalog_rho_x_prop from last loop, saves calls to rho
                !computing average acceptance rate, will want to have way of outputing at some point
                a_ave = (1.0_dp/loop_index)*(a_ave*(loop_index-1.0_dp)+1.0_dp)
            else
                a_ave = (1.0_dp/loop_index)*(a_ave*(loop_index-1))
            end if
            if (loop_index.ge.n_burned) then
                thinning_counter = thinning_counter+1
            end if
            if (thinning_counter .ge. thinning_interval) then
                out_array_index = out_array_index+1
                if (out_array_index>(n_steps-n_burned)/thinning_interval) then
                    e_code = e_code+2**2
                    print*, 'ERROR: Output larger than output array', loop_index
                    exit
                end if
                if (a.ne.a) then
                    e_code=e_code+1
                    print*, 'ERROR: NAN acceptance at:', loop_index
                    EXIT
                end if
                if (a_ave.le.0) then
                    if (a_ave<0) then
                        e_code=e_code+2**1
                        print*, 'ERROR: Acceptance rate negative at:', loop_index
                        EXIT
                    end if
                    print*, 'WARNING: Acceptance rate 0 at:', loop_index
                end if
                samples(out_array_index,:) = x
                thinning_counter=0
            end if
        end do


        if (present(average_accept)) then
            average_accept = a_ave
        end if

    end subroutine mcmc_sample

    !> @brief Runs a version of mcmc_sample but every adapt_interval steps adjusts s, based on an exponetial average of the accpetance rate
    subroutine mcmc_adapt(s_out, log_rho, x_0, n_steps, s_0, e_code,&
         s_max, s_min, memory, adapt_interval, dof_coefficients, density_dimension, seed)
        implicit none
        real(dp), intent(out) :: s_out !>adapted step size at end of run
        procedure(log_rho_interface) :: log_rho !>the (log) density to be sampled from, takes input, x, and parametrisation dof_coefficients
        integer, intent(in) :: density_dimension !>size of x
        integer, intent(in) :: n_steps !>total number of steps to take
        integer, intent(in) :: adapt_interval !>interval between performing the step size adaption routine
        integer, intent(in), dimension(:), optional :: seed !>the random seed, will use current state if not set
        real(dp), dimension(:), intent(in) :: dof_coefficients!> the parameterisation of log_rho
        real(dp), intent(in) :: s_0 !>intitial guess at step size
        real(dp), intent(in) :: s_max, s_min !>bounds for search space in s
        real(dp), intent(in) :: memory !>average memory of exponetial kernal
        real(dp), dimension(density_dimension), intent(in) :: x_0 !>a starting point for the chain
        integer, intent(out) :: e_code !>an error code, 0 is successs, 1, is a NaN accpetance, 2 is a negative acceptance, 4 is array size mismatch
        real(dp), dimension(density_dimension) :: x,x_prop
        real(dp) :: a,a_unnorm,norm_cnst,a_ave,rand_sample,log_rho_x,log_rho_x_prop,s, factor, a_ave_
        real(dp) :: a_min, a_max !keep in close range about 0.235
        integer :: loop_index, adapt_count

        a_min = 0.234_dp
        a_max = 0.236_dp
        x = x_0
        a_unnorm = 0.0_dp
        log_rho_x = log_rho(x_0,dof_coefficients)
        s=s_0
        norm_cnst=1.0_dp
        a_ave_=1.0_dp
        factor=sqrt(2.0_dp)

        !uses current random state if none is given
        if (present(seed)) then
            call random_seed(put=seed)
        end if

        do loop_index=1,n_steps
            x_prop = mcmc_propose(x,s,density_dimension)
            call mcmc_accept(a, x, x_prop, log_rho, s, log_rho_x, log_rho_x_prop,dof_coefficients)
            a = exp(a)
            call random_number(rand_sample)
            if (a>rand_sample) then
                x = x_prop
                log_rho_x = log_rho_x_prop !using datalog_rho_x_prop from last loop, saves calls to rho
                !computing average acceptance rate, will want to have way of outputing at some point
                a_unnorm = 1.0_dp + exp(-1.0_dp/memory)*a_unnorm
                a_ave_= (1.0_dp+(loop_index-1)*a_ave_)/loop_index
            else
                a_unnorm = exp(-1.0_dp/memory)*a_unnorm
                a_ave_= ((loop_index-1)*a_ave_)/loop_index
            end if
            norm_cnst = 1.0_dp + exp(-1.0_dp/memory)*norm_cnst

            adapt_count = adapt_count+1
            if (adapt_count .ge. adapt_interval) then
                if (a_ave_.le.0) then
                    if (a_ave_<0) then
                        e_code=e_code+2**1
                        print*, 'ERROR: Acceptance rate negative at:', loop_index
                        EXIT
                    end if
                    print*, 'WARNING: Acceptance rate 0 at:', loop_index
                end if
                adapt_count=0
                a_ave = a_unnorm/norm_cnst
                if (a_ave<a_min) then
                    s=s/factor
                    if (s<s_min) then
                        s=s_min
                    end if
                else if (a_ave>a_max) then
                    s=s*factor
                    if (s>s_max) then
                        s=s_max
                    end if
                end if
            end if
        end do

        s_out=s
        !e_code 0 => no error
        e_code=0

    end subroutine mcmc_adapt

  end module
