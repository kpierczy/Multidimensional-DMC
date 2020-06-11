%========================================================================
% File name   : DMC.m
% Date        : 28th April 2020
% Version     : 1.0.0
% Author      : Krzysztof Pierczyk
% Description : Generic class implementing DMC MIMO regulator. Class
%               takes into account both controller (CVs) and uncontroller
%               (DVs) inputs.
%
%               Regulator's construction is based on the given DMC's
%               parameters and handle to a function simulating
%               a regulated object. Step response of all tracks are
%               gathered automatically under controll of a few
%               parameters in the structure passed as argument to class'
%               constructor.
%
%               Regulator's parameters as well as regulated object can
%               be modified 'in the fly'.
%
%               Regulator can be configured to take into account 
%               limitations of input signals.
%
% To do list  :
%
%    (1) Introduce internall nu and ny values instead of using
%        size(obj.step_responses, ...)
%    (2) Introduce 'fi' parameters matrix and extend obj.lambda
%        to matrxi of parameters
% 
%========================================================================

classdef DMC < handle
    
    properties (Access = public, Dependent)
        
        % DMC regulator's parameters (interface)
        N      (1, 1) uint32 {mustBeInteger, mustBePositive}
        Nu     (1, 1) uint32 {mustBeInteger, mustBePositive}
        D      (1, 1) uint32 {mustBeInteger, mustBePositive}
        Dz     (1, 1) uint32 {mustBeInteger, mustBeNonnegative}
        lambda (1, 1) double {mustBeReal   , mustBeNonnegative}
        
    end
    
    properties(SetAccess = private, GetAccess = public)
       
        % DMC matrices
        K   (:, :) double {mustBeReal}
        Mp  (:, :) double {mustBeReal}
        Mzp (:, :) double {mustBeReal}
        
    end
    
    properties (Access = public)
        
        % Limits of controlled variables
        umin  (:, 1) double {mustBeReal} = -Inf
        umax  (:, 1) double {mustBeReal} =  Inf
        dumin (:, 1) double {mustBeReal} = -Inf
        dumax (:, 1) double {mustBeReal} =  Inf

        
        
        % @note : DMC is an incremental regulator (it computes CVs
        %         increases, not absolute values). To minimize number
        %         of arguments passed to the 'compute()' method every
        %         iteration of controll loop DMC class implements
        %         internal memory region to keep track of regulation
        %         history.
        %
        %         Following properties are public for cases when 
        %         regulation process is switched to manual for some
        %         period of time. When regulation is sitched-back to
        %         automatic regulation these parameters should be updated
        %         to ensure smooth start-up.
        
        % Last CVs values
        uk_1  (:, 1) double {mustBeReal}
        
        % Last DVs values
        zk_1  (:, 1) double {mustBeReal}
        
        % Vector of differences of subsequent CVs on the steerage horizon.
        dUp   (:, 1) double {mustBeReal} 
        
        % Vector of differences of subsequent DVs on the steerage horizon.
        dZp  (:, 1) double {mustBeReal} 
        
    end
    
    % Internal interface elements
    properties (Access = private)
        
        % Initialization flag
        initialized = false
        
        % DMC regulator's parameters (internal copies)
        N_shadow      (1, 1) uint32 {mustBeInteger, mustBePositive}    = 1
        Nu_shadow     (1, 1) uint32 {mustBeInteger, mustBePositive}    = 1
        D_shadow      (1, 1) uint32 {mustBeInteger, mustBePositive}    = 1
        Dz_shadow     (1, 1) uint32 {mustBeInteger, mustBeNonnegative} = 0
        lambda_shadow (1, 1) double {mustBeReal   , mustBeNonnegative} = 0
        
        % DMC matrices (internal copies)
        K_shadow   (:, :) double {mustBeReal}
        Mp_shadow  (:, :) double {mustBeReal}
        Mzp_shadow (:, :) double {mustBeReal}
        
        % Recently saved step responses
        step_responses (:, :, :) double {mustBeReal}
        
        % Recently saved step responses of the output-disturbance tracks
        z_step_responses (:, :, :) double {mustBeReal}
        
    end
    
    
    %====================================================================
    %%%%%%%%%%%%%%%%%%%%%%%%%% Public methods %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %====================================================================
    
    methods (Access = public)
        
        %================================================================
        % DMC regulator initialization. 
        % 
        % @param step_responses_struct : Structure containing information
        %        used to gather all required step responses. 
        %        The structure should containt following field:
        %
        %        .object - handle to the "function Y = f(U)"
        %                  or "function Y = f(U, Z) function representing 
        %                  regulated object. It should take two 2D arrays
        %                  representing CVs and DVs signals (matrices'
        %                  rows refer to a single CV/DV). It should return 
        %                  output for moments from k = 1 to
        %                  k = (size(input))
        %
        %                  If number of columns of U and Z matrices are
        %                  not equal, function should throw an error.
        %
        %        .shape  - [nu, nz, ny] array containing object's 
        %                  dimensions
        %
        %        .init_point_u - vector of (nu x 1) shape conatining
        %                  the initial CVs work point used to perform
        %                  steps during step responses' gathering
        %
        %        .init_point_z - vector of (nz x 1) shape conatining
        %                  the initial DVs work point used to perform
        %                  steps during step responses' gathering
        %
        %        .step_size_u - vector of (nu x 1) shape conatining
        %                  steps sizes of CVs used to gather step
        %                  responses
        %
        %        .step_size_z - vector of (nz x 1) shape conatining
        %                  steps sizes of DVs used to gather step
        %                  responses
        %
        %        .tol (optional, defaul = 0.0001) - maximum difference
        %                  between two values (e.g. subsequent output
        %                  values) that can be approximated as 0
        %
        %        .sim_time (optional, default = 1000) - length of the
        %                  simulation used to get step responses
        %
        % @note : Function will thrown an error for all objects 
        %         that stabilize longer than 'sim_time'. An object 
        %         is considered stable if absolute differences between 
        %         two subsequent values of all outputs are smaller than 
        %         .tol
        %
        %
        % @param param_struct : structure containing DMC parameters
        %
        %        .N - prediction horizon; positive integer value.
        %             If .N is greater than dynamic rank .D (see below)
        %             it is set to .D
        %
        %        .Nu - control horizon; positive integer value.
        %             If .Nu is greater than dynamic rank .D (see below)
        %             it is set to .D
        %
        %        .D - dynamic rank; positive integer value.
        %             If .D is greater than stabilization period of
        %             the slowest CV-PV track it is cut to the value
        %             of this period
        %
        %        .Dz - disturbance dynamic rank; non-negative integer 
        %             value. If .Dz is greater than stabilization period 
        %             of the slowest DV-PV track it is cut to the value
        %             of this period.
        %             If .Dz is set to 0 DV-PV tracks are not taken
        %             into account in regulation process
        %
        %         .lambda - non-negative real value controlling 
        %             punishment for the CVs changes during optimisation
        %             process
        %================================================================
        function obj = DMC(step_responses_struct, param_struct)
            
            % Initial parameters verification
            mustBeInteger(param_struct.N);
            mustBePositive(param_struct.N);
            mustBeInteger(param_struct.Nu);
            mustBePositive(param_struct.Nu);
            mustBeInteger(param_struct.D);
            mustBePositive(param_struct.D);
            mustBeInteger(param_struct.Dz);
            mustBeNonnegative(param_struct.D);
            mustBeReal(param_struct.lambda);
            mustBeNonnegative(param_struct.lambda);
            
            
            % Gather step responses
            obj.update_model(step_responses_struct)
            
            % Initialize DMC parameters
            obj.D_shadow      = param_struct.D; 
            obj.Dz_shadow     = param_struct.Dz;
            obj.N_shadow      = param_struct.N;
            obj.Nu_shadow     = param_struct.Nu;
            obj.lambda_shadow = param_struct.lambda;
            obj.verify_params();
            
            % Reset limits and memory after parameters asignment
            obj.reset();
            
            % Mark object initialized
            obj.initialized = true;
            
            % Compute DMC matrices
            obj.solve();
            
        end
        
        %================================================================
        % Function changes saved step responses basing on the new model 
        % of the object.
        %
        % @param step_responses_struct : @see DMC(...)
        %
        % @note : If DMC parameters doesn't meet requirements described
        %         in the constructor's commet they are modified as
        %         described above (@see DMC(...))
        %
        % @note : Regulator's memory (i.e .uk_1, .dUp and .dZp vectors)
        %         are set to zero and should be updated after method's
        %         call if required. All limits are discarded.
        %================================================================
        function update_model(obj, step_responses_struct)
            
            % Alias argument's name
            srs = step_responses_struct;
            
            % Compute step responses for CV-PV tracks
            step_responses_t   = cell(srs.shape(1), srs.shape(3));
            dynamic_rank = 0;
            for u = 1:srs.shape(1)
                for y = 1:srs.shape(3)
                    
                    % Prepare SISO bind
                    bind_struct.object     = srs.object;
                    bind_struct.shape      = srs.shape;
                    bind_struct.U          = srs.init_point_u;
                    if isfield(srs, 'init_point_z')
                        bind_struct.Z      = srs.init_point_z;
                    else
                        bind_struct.Z      = [];
                    end                    
                    bind_struct.input_num  = u;
                    bind_struct.output_num = y;
                    bind_struct.input_type = 'CV';
                    
                    % Initialize structure for single step response calculation
                    srs_SISO.object = @(input)(DMC.bind(bind_struct, input));
                    srs_SISO.init_point = srs.init_point_u(u);
                    srs_SISO.step_size  = srs.step_size_u(u);
                    if isfield(srs, 'tol')
                        srs_SISO.tol = srs.tol;
                    end
                    if isfield(srs, 'sim_time')
                        srs_SISO.sim_time = srs.sim_time;
                    end
                    
                    % Compute step response
                    step_responses_t{u, y} = DMC.step_response_SISO(srs_SISO);
                    
                    % Update dynamic_rank
                    if size(step_responses_t{u, y}, 2) > dynamic_rank
                        dynamic_rank = size(step_responses_t{u, y}, 2);
                    end
                    
                end
            end
            
            % Extend step responses to a dynamic rank value
            obj.D_shadow = dynamic_rank;
            obj.step_responses = zeros(srs.shape(1), srs.shape(3), obj.D);
            for u = 1:srs.shape(1)
                for y = 1:srs.shape(3)
                    length = size(step_responses_t{u, y}, 2);
                    obj.step_responses(u, y, 1:length) = step_responses_t{u, y};
                    obj.step_responses(u, y, length:end) = obj.step_responses(u, y, length);
                end
            end
            
            % Check if DVs are taken into account
            if srs.shape(2) ~= 0

                % Compute step responses for DV-PV tracks
                z_step_responses_t = cell(srs.shape(2), srs.shape(3));
                z_dynamic_rank = 0;
                for z = 1:srs.shape(2)
                    for y = 1:srs.shape(3)
                        
                        % Prepare SISO bind
                        bind_struct.object     = srs.object;
                        bind_struct.shape      = srs.shape;
                        bind_struct.U          = srs.init_point_u;
                        bind_struct.Z          = srs.init_point_z;
                        bind_struct.input_num  = z;
                        bind_struct.output_num = y;
                        bind_struct.input_type = 'DV';
                        
                        % Initialize structure for single step response calculation
                        srs_SISO.object = @(input)(DMC.bind(bind_struct, input));
                        srs_SISO.init_point = srs.init_point_z(z);
                        srs_SISO.step_size  = srs.step_size_z(z);
                        if isfield(srs, 'tol')
                            srs_SISO.tol = srs.tol;
                        end
                        if isfield(srs, 'sim_time')
                            srs_SISO.sim_time = srs.sim_time;
                        end
                        
                        % Compute step response
                        z_step_responses_t{z, y} = DMC.step_response_SISO(srs_SISO);
                        
                        % Update dynamic_rank
                        if size(z_step_responses_t{z, y}, 2) > z_dynamic_rank
                            z_dynamic_rank = size(z_step_responses_t{z, y}, 2);
                        end
                        
                    end
                end
                
                % Extend step responses to a disturbance dynamic rank value
                obj.Dz_shadow = z_dynamic_rank;
                obj.z_step_responses = zeros(srs.shape(2), srs.shape(3), obj.Dz);
                for z = 1:srs.shape(2)
                    for y = 1:srs.shape(3)
                        length = size(z_step_responses_t{z, y}, 2);
                        obj.z_step_responses(z, y, 1:length) = z_step_responses_t{z, y};
                        obj.z_step_responses(z, y, length:end) = obj.z_step_responses(z, y, length);
                    end
                end

            else
                
                obj.Dz_shadow = 0;
                obj.z_step_responses = [];

            end
            
            % Reset some parameters if object was used yet
            if obj.initialized
                obj.verify_params();
                obj.solve();
                obj.reset();
            end
            
        end
        
        %================================================================
        % Computes optimal CV value for the given process' output and
        % desired output on the prediction horizon.
        %
        % @param y_zad : matrix of size (ny x X) representing desired 
        %        PVs values for next X discrete moments. Rows of the
        %        matrix refer to subsequent outputs when columns mark
        %        discrete moments.
        %        If X < N the Y_zad vector is extended to the
        %        X * ny length with the last element of the trajectory
        %        for each PV.
        %        If X > N the Y_zad vector is cut to X*ny.
        %
        % @param y : vector of actual PVs values
        %
        % @param z : vector of DVs
        %================================================================
        function u = compute(obj, y_zad, y, z)
            
            % Check vector sized
            if size(y_zad, 1) ~= size(obj.step_responses, 2)
               error('y_zad vector is of the wrong size'); 
            end
            if size(y, 1) ~= size(obj.step_responses, 2)
               error('y vector is of the wrong size'); 
            end
            
            % If noise is taken into account update dZp property
            if obj.Dz ~= 0 
                if exist('z', 'var')

                    % Verify z vector size
                    if size(z, 1) ~= size(obj.z_step_responses, 1)
                       error('z vector is of the wrong size'); 
                    end

                    % Circshift dZp vector
                    obj.dZp = circshift(obj.dZp, size(obj.z_step_responses, 1));

                    % Compute DVs  inreases
                    obj.dZp(1:size(obj.z_step_responses, 1)) = z - obj.zk_1;
                    
                end
            end
                
            % Initialize reshaped y_zad vector
            Y_zad = zeros(obj.N * size(obj.step_responses, 2), 1);
            
            % Reshape and cut y_zad vector
            if size(y_zad, 2) > obj.N
                
                % Fill reshaped y_zad vector
                for i = 1:obj.N
                   Y_zad( ...
                       (i - 1) * size(obj.step_responses, 2) + 1 : ...
                            i  * size(obj.step_responses, 2)       ...
                    ) = y_zad(:, i);
                       
                end

            % Reshape y_zad with or without extension
            else
                
                % Reshape y_zad
                for i = 1:size(y_zad, 2)
                   Y_zad( ...
                       (i - 1) * size(obj.step_responses, 2) + 1 : ...
                            i  * size(obj.step_responses, 2)       ...
                    ) = y_zad(:, i);
                end
                
                % Extend y_zad
                if size(y_zad, 2) < obj.N
                    for i = size(y_zad, 2):obj.N
                        Y_zad( ...
                            (i - 1) * size(obj.step_responses, 2) + 1 : ...
                                 i  * size(obj.step_responses, 2)       ...
                         ) = y_zad(:, end);
                    end
                end
            end
            
            % Validate y vector
            if size(y, 1) ~= size(obj.step_responses, 2)
               error('Too few output values given!') 
            end
            
            % Create Y vector
            Y = repmat(y, obj.N, 1);
        
            %---------------- Regulation law computation -----------------
            
            % Compute unrestrained trajectory
            if obj.Dz > 0
                Y_0 = Y + obj.Mp * obj.dUp + obj.Mzp * obj.dZp;
            else
                Y_0 = Y + obj.Mp * obj.dUp;
            end
            
            % Compute CVs increases
            dU = obj.K(1:size(obj.step_responses, 1), :) * (Y_zad - Y_0);
            
            %------------------------------------------------------------
            
            % Get CVs increases for this iteration only
            du = dU(1:size(obj.step_responses, 1));            
            
            % Apply dU value constraints
            for i = 1:size(du)
                if du(i) > obj.dumax(i)
                    du(i) = obj.dumax(i);
                elseif du(i) < obj.dumin(i)
                    du(i) = obj.dumin(i);
                end
            end

            % Circshift dUp vector
            obj.dUp = circshift(obj.dUp, size(obj.step_responses, 1));

            % Save CVs  inreases
            obj.dUp(1:size(obj.step_responses, 1)) = du;
            
            % Compute CVs
            u = obj.uk_1 + du;
            
            % Apply U value constraints
            for i = 1:size(obj.step_responses, 1)
                if u(i) > obj.umax(i)
                    u(i) = obj.umax(i);
                elseif u(i) < obj.umin(i)
                    u(i) = obj.umin(i);
                end
            end
            
            % Save CVs and DVs values for next iteration
            obj.uk_1 = u;
            if exist('z', 'var')
                obj.zk_1 = z;
            end
            
        end    
        
        
        %================================================================
        % Resets regulator's limits and memory to the default state.
        %================================================================
        function reset(obj)
           
            % Clear memory
            obj.uk_1 = zeros(size(obj.step_responses  , 1), 1);
            obj.zk_1 = zeros(size(obj.z_step_responses, 1), 1);
            obj.dUp  = zeros(size(obj.step_responses  , 1) * (obj.D - 1), 1);
            obj.dZp  = zeros(size(obj.z_step_responses, 1) * obj.Dz     , 1);

            % Reset limits of controlled variables
            obj.umin   = ones(size(obj.step_responses, 1), 1) * (-Inf);
            obj.umax   = ones(size(obj.step_responses, 1), 1) * Inf;
            obj.dumin  = ones(size(obj.step_responses, 1), 1) * (-Inf);
            obj.dumax  = ones(size(obj.step_responses, 1), 1) * Inf;
            
        end
        
    end
    
    
    %====================================================================
    %%%%%%%%%%%%%%%%%%%%%%%%%% Private methods %%%%%%%%%%%%%%%%%%%%%%%%%%
    %====================================================================
    
    methods (Access = private)
        
        %================================================================
        % Verifies DMC's parameters after updating the model of the 
        % object.
        %
        % @returns : True if some parameters were modified
        %================================================================
        function modified = verify_params(obj)
            
            % Modification flag
            modified = false;
                        
            % Initialize DMC parameters
            if obj.D > size(obj.step_responses, 3)
                obj.D_shadow =  size(obj.step_responses, 3);
                modified = true;
            end
            
            if obj.Dz > size(obj.z_step_responses, 3)
                obj.Dz_shadow = size(obj.z_step_responses, 3);
                modified = true;
            elseif obj.Dz > 0 && size(obj.z_step_responses, 1) == 0
                obj.Dz_shadow = 0;
                modified = true;
            end
            
            if obj.N > obj.D
                obj.N_shadow =  obj.D;
                modified = true;
            end
            
            if obj.Nu > obj.D
                obj.Nu_shadow =  obj.D;
                modified = true;
            end        
            
        end
        
        %================================================================
        % Computes analitical solution of the DMC optimization problem 
        % returning K, Mp and Mzp matrices used to calculate DMC 
        % algorithm's output.
        %================================================================
        function solve(obj)

            % Initialize M matrix
            M  = zeros( ...
                    size(obj.step_responses, 2) * obj.N, ...
                    size(obj.step_responses, 1) * obj.Nu ...
                 );
            for i = 1:obj.N
                for j = 1:obj.Nu
                    if i >= j
                        
                        % Create S matrix
                        S = obj.S(i - j + 1);
                        
                        % Fill field in M matrix
                        M((i-1)*size(obj.step_responses, 2) + 1 : i*size(obj.step_responses, 2), ...
                          (j-1)*size(obj.step_responses, 1) + 1 : j*size(obj.step_responses, 1)) = S;
                      
                    end
                end
            end

            % Initialize Mp matrix
            obj.Mp = zeros( ...
                        size(obj.step_responses, 2) * obj.N,     ...
                        size(obj.step_responses, 1) *(obj.D - 1) ...
                     );
            for j = 1:(obj.D - 1)
                
                % Create Sj matrix
                Sj = obj.S(j);
                
                for i = 1:obj.N
                    
                    % Create S matrix
                    S = obj.S(i+j);
                    
                    % Fill Mp's field with Ss' difference
                    obj.Mp((i-1)*size(obj.step_responses, 2) + 1 : i*size(obj.step_responses, 2), ...
                           (j-1)*size(obj.step_responses, 1) + 1 : j*size(obj.step_responses, 1)) = S - Sj;

                end
            end
            

            % Initialize Mzp matrix
            if obj.Dz ~= 0
                obj.Mzp = zeros(....
                              size(obj.step_responses,   2) * obj.N, ...
                              size(obj.z_step_responses, 1) * obj.Dz ...
                          );
               for j = 1:obj.Dz

                    % Create Sj matrix
                    if j > 1
                        Sj = obj.Sz(j-1);
                    end

                    for i = 1:obj.N

                        % Create S matrix
                        S = obj.Sz(i+j-1);

                        % Fill Mzp's field with Ss' difference
                        if j > 1
                            obj.Mzp((i-1)*size(obj.step_responses, 2) + 1 : i*size(obj.step_responses, 2), ...
                                    (j-1)*size(obj.step_responses, 1) + 1 : j*size(obj.step_responses, 1)) = S - Sj;
                        else
                            obj.Mzp((i-1)*size(obj.step_responses, 2) + 1 : i*size(obj.step_responses, 2), ...
                                    (j-1)*size(obj.step_responses, 1) + 1 : j*size(obj.step_responses, 1)) = S;
                        end

                    end
               end
            else
                obj.Mzp = [];
            end
             
            % Solve optimization problem
            obj.K = (M' * M + obj.lambda * eye(size(obj.step_responses, 1) * obj.Nu))^(-1) * M';
        end
        
        %================================================================
        % Auxiliary function that constructs S_l matrix based on the
        % l index
        %
        % @param l : index of the samples
        %================================================================
        function S_l = S(obj, l)
            S_l = zeros(size(obj.step_responses, 2), size(obj.step_responses, 1));
            if l <= obj.D
                for y = 1:size(obj.step_responses, 2)
                    for u = 1:size(obj.step_responses, 1)
                        S_l(y, u) = obj.step_responses(u, y, l);
                    end
                end
            else
                for y = 1:size(obj.step_responses, 2)
                    for u = 1:size(obj.step_responses, 1)
                        S_l(y, u) = obj.step_responses(u, y, obj.D);
                    end
                end
            end
        end
        
        %================================================================
        % Auxiliary function that constructs Sz_l matrix based on the
        % l index
        %
        % @param l : index of the samples
        %================================================================
        function Sz_l = Sz(obj, l)
            Sz_l = zeros(size(obj.z_step_responses, 2), size(obj.z_step_responses, 1));
            if l <= obj.Dz
                for y = 1:size(obj.z_step_responses, 2)
                    for z = 1:size(obj.z_step_responses, 1)
                        Sz_l(y, z) = obj.z_step_responses(z, y, l);
                    end
                end
            else
                for y = 1:size(obj.z_step_responses, 2)
                    for z = 1:size(obj.z_step_responses, 1)
                        Sz_l(y, z) = obj.z_step_responses(z, y, obj.Dz);
                    end
                end
            end
        end
        
    end
    
    methods (Access = private, Static = true)
        
        %================================================================
        % Computes step response suitable to use with DMC algorithm
        % for a single input - single output track.
        %
        % @param step_responses_struct : Structure containing information
        %        used to gather all required step responses. 
        %        The structure should containt following field:
        %
        %        .object - handle to the "function output = object(input)" 
        %                  function representing regulated object.
        %                  It should take an input vector representing
        %                  input values during simulation time and should 
        %                  return output for moments from k = 1 to 
        %                  k = (size(input))
        %
        %        .init_point - initial input value used to perform
        %                  step during gathering of a step response 
        %
        %        .step_size - input step size used to perform gather
        %                  a step response 
        %
        %        .tol (optional, defaul = 0.0001) - maximum difference
        %                  between two values (e.g. subsequent output
        %                  values) that can be approximated as 0
        %
        %        .sim_time (optional, default = 1000) - length of the
        %                  simulation used to get step responses
        %          
        % @note : Function will thrown an error for all objects 
        %         that stabilize longer than 'sim_time'. An object 
        %         is considered stable if absolute difference between 
        %         two subsequent values the output is smaller than 
        %         .tol
        %================================================================
        function Y = step_response_SISO(step_response_struct)

            % Alias argument's name
            srs = step_response_struct;
            
            % Check optional arguments
            if isfield(srs, 'tol')
                mustBePositive(srs.tol);
                tol = srs.tol;
            else
                tol = 0.0001;
            end
            
            if isfield(srs, 'sim_time')
                mustBeInteger(srs.sim_time);
                mustBePositive(srs.sim_time);
                sim_time = srs.sim_time;
            else
                sim_time = 1000;
            end


            % Search the work point
            U = ones(1, sim_time) * srs.init_point;
            Y = srs.object(U);

            % Get a stabilization time and work point
            stabilization_time = find(abs(Y - Y(end)) < tol, 1, 'first');
            if isempty(stabilization_time)
               error("Object's track could not be stabilized!") 
            end
            Ypp = Y(stabilization_time);

            % Measure step response
            U(stabilization_time:end) = step_response_struct.init_point + step_response_struct.step_size;
            Y = step_response_struct.object(U);

            % Scale step response
            Y = (Y - Ypp) / step_response_struct.step_size;
            
            % Trim step response from left and right
            if stabilization_time ~= size(Y, 2)
                Y = Y(stabilization_time+1 : end);
            else
                error('Object could not be stabilized!')
            end
            Y = Y(1 : find(abs(Y - Y(end)) < tol, 1, 'first'));
            
        end

        %================================================================
        % Binds MIMO object function to a single SISO track with given
        % parameters.
        %
        % @param bind_struct : structure that describes way bind should
        %        be performed. Structure's fields are:
        %
        %        .object - binded MIMO function in a "function Y = f(U)"
        %                  or "function Y = f(U, Z)" form. To read about
        %                  details @see DMC(...)
        %
        %        .shape  - [nu, nz, ny] array containing object's 
        %                  dimensions
        %
        %        .U - vector of 'nu' elements conatining values of
        %             CVs that will be bound
        %
        %        .Z - vector of 'nz' elements conatining values of
        %             DVs that will be bound        
        %
        %        .input_num - number of the free input
        %
        %        .output_num - number of the required PV
        %
        %        .input_type - type of the free input ('CV' or 'DV') 
        %
        % @param input : free input course
        %================================================================
        function output = bind(bind_struct, input)
            
            % Initialize input
            U = zeros(bind_struct.shape(1), size(input, 2));
            Z = zeros(bind_struct.shape(2), size(input, 2));
            
            % Fill input matrices
            for u = 1:size(U, 1)
                U(u, :) = bind_struct.U(u);
            end
            for z = 1:size(Z, 1)
                Z(z, :) = bind_struct.Z(z);
            end
            
            % Set free CV
            if bind_struct.input_type == 'CV'
                U(bind_struct.input_num, :) = input;
            elseif bind_struct.input_type == 'DV'
                Z(bind_struct.input_num, :) = input;
            else
                error('Wrong input type!')
            end
            
            % Compute bind
            if bind_struct.shape(2) ~= 0
                output = bind_struct.object(U, Z);
            else
                output = bind_struct.object(U);
            end
            output = output(bind_struct.output_num, :);
            
        end
        
    end
    
    
    %====================================================================
    %%%%%%%%%%%%%%%%%%%%%%%% Getters & setters %%%%%%%%%%%%%%%%%%%%%%%%%%
    %====================================================================
    
    methods
       
        %================================================================
        % Block of getters & setters used to remove the need of explicit
        % DMC matrices calculation.
        %================================================================
        
        function val = get.N(obj)
            val = obj.N_shadow;
        end
        function set.N(obj, N)
            obj.N_shadow = N;
            obj.verify_params()
            obj.solve();
            obj.reset();
        end
        
        function val = get.Nu(obj)
            val = obj.Nu_shadow;
        end
        function set.Nu(obj, Nu)
            obj.Nu_shadow = Nu;
            obj.verify_params()
            obj.solve();
            obj.reset();
        end
        
        function val = get.D(obj)
            val = obj.D_shadow;
        end
        function set.D(obj, D)
            obj.D_shadow = D;
            obj.verify_params()
            obj.solve();
            obj.reset();
        end
        
        function val = get.Dz(obj)
            val = obj.Dz_shadow;
        end
        function set.Dz(obj, Dz)
            obj.Dz_shadow = Dz;
            obj.verify_params()
            obj.solve();
            obj.reset();
        end
        
        function val = get.lambda(obj)
            val = obj.lambda_shadow;
        end
        function set.lambda(obj, lambda)
            obj.lambda_shadow = lambda;
            obj.verify_params()
            obj.solve();
            obj.reset();
        end
        
        function K = get.K(obj)
            K = obj.K_shadow;
        end
        function set.K(obj, K)
            obj.K_shadow = K;
        end        
        
        function Mp = get.Mp(obj)
            Mp = obj.Mp_shadow;
        end
        function set.Mp(obj, Mp)
            obj.Mp_shadow = Mp;
        end
        
        function Mzp = get.Mzp(obj)
            Mzp = obj.Mzp_shadow;
        end
        function set.Mzp(obj, Mzp)
            obj.Mzp_shadow = Mzp;
        end
        
    end
end

