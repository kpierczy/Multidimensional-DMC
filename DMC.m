%========================================================================
% File name   : DMC.m
% Date        : 11th June 2020
% Version     : 1.2.0
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
%========================================================================

classdef DMC < handle
    
    properties (Access = public, Dependent)
        
        % DMC regulator's parameters (interface)
        N      (1, 1) uint32 {mustBeInteger, mustBePositive}
        Nu     (1, 1) uint32 {mustBeInteger, mustBePositive}
        D      (1, 1) uint32 {mustBeInteger, mustBePositive}
        Dz     (1, 1) uint32 {mustBeInteger, mustBeNonnegative}
        psi    (:, :) double {mustBeReal   , mustBeNonnegative}
        lambda (:, :) double {mustBeReal   , mustBeNonnegative}
        
        % note : 'psi' and 'lambda' can be either scalar or diagonal
        %        matrices. If they are matrices they have to be of the
        %        (N*ny) x (N*ny) and (Nu*nu) x (Nu*nu) shape respectively.
        %        Otherwise, they are extended to they eye(N*ny) and eye(Nu*nu)
        %        matrices before solving DMC optimisation problem
        
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
    
    properties (SetAccess = private, GetAccess = public, Dependent)
        
        % DMC matrices
        K   (:, :) double {mustBeReal}
        Mp  (:, :) double {mustBeReal}
        Mzp (:, :) double {mustBeReal}
        
        % View on the internal step responses
        step_responses   (:, :, :) double {mustBeReal}
        z_step_responses (:, :, :) double {mustBeReal}
        
    end
    
    properties(SetAccess = private, GetAccess = public)
               
        % Object's shape
        nu uint32 {mustBeNonnegative} = 0
        nz uint32 {mustBeNonnegative} = 0
        ny uint32 {mustBeNonnegative} = 0
        
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
        lambda_shadow (:, :) double {mustBeReal   , mustBeNonnegative} = 0
        psi_shadow    (:, :) double {mustBeReal   , mustBeNonnegative} = 0
        
        % DMC matrices (internal copies)
        K_shadow   (:, :) double {mustBeReal} = 0
        Mp_shadow  (:, :) double {mustBeReal} = 0
        Mzp_shadow (:, :) double {mustBeReal} = 0
        
        % Recently saved step responses
        step_responses_shadow (:, :, :) double {mustBeReal} = 0
        
        % Recently saved step responses of the output-disturbance tracks
        z_step_responses_shadow (:, :, :) double {mustBeReal} = 0
        
    end
    
    
    %====================================================================
    %%%%%%%%%%%%%%%%%%%%%%%%%% Public methods %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %====================================================================
    
    methods (Access = public)
        
        %================================================================
        % DMC regulator initialization. 
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
        % @note : calculating optimal DMC matrices requires 'lambda' and
        %         'psi' diagonal matrices containing weights for adjusting
        %         balance between importance of particular output's
        %         stabilization and variance of the CVs. These are by
        %         default set to eye(Nu * nu) and eye(N * ny). It is
        %         because given horizonts (N, Nu) can turn out to be greater
        %         that the maximal time of outputs stabilization during
        %         model's initialization process. In this case they are
        %         trimmed so given 'lambda' and 'psi' matrices would be
        %         of the wrong size. 
        %         After initialization of the internal model actual N
        %         and Nu values can be checked and proper materices can
        %         be set using '.' operator
        %
        % @note : parameters are set in constructor but object remains
        %         uninitialized up to the first call of updateModel() or
        %         updateModelFunction() method
        %
        % @note : class does not controll size of the .lambda and .psi
        %         matrices. They should be set to appropriate sizes
        %         by the user
        %================================================================
        function obj = DMC(param_struct)
            
            % Initial parameters verification
            mustBeInteger(param_struct.N);
            mustBePositive(param_struct.N);
            mustBeInteger(param_struct.Nu);
            mustBePositive(param_struct.Nu);
            mustBeInteger(param_struct.D);
            mustBePositive(param_struct.D);
            mustBeInteger(param_struct.Dz);
            mustBeNonnegative(param_struct.D);
            
            % Initialize DMC parameters
            obj.D_shadow      = param_struct.D; 
            obj.Dz_shadow     = param_struct.Dz;
            obj.N_shadow      = param_struct.N;
            obj.Nu_shadow     = param_struct.Nu;
            obj.lambda_shadow = 1;
            obj.psi_shadow    = 1;
            
        end
        
        %================================================================
        % Function (re)initializes regulator using set of the step
        % responses. Responses can be raw signals taken from the object.
        % If so, they are scaled by the function wth respect to given
        % scaling parameters).
        %
        % @param stepResponsesStruct - structure containing all
        %        informations eqruired to initialize DMC step responses:
        %
        %        .stepResponses - cell aray of size (nu x ny). Each cell
        %             contains step response of the corresponding input-
        %             output track. A step response is the column vector.
        %
        %        .zStepResonses - cell aray of size (nz x ny). Each cell
        %             contains step response of the corresponding 
        %             disturbance - output track. A step response is the
        %             column vector.
        %
        %        .workPoint (optional, default 0) - vector of size
        %             (ny x 1) containng values of the outputs in the work
        %             point that step responses were acquired from
        %
        %         .stepSize (optional, default 1) - vector of size
        %             (nu x 1) containng values of the inputs steps used
        %             to gather step responses
        %
        %         .zStepSize (optional, default 1) - vector of size
        %             (nz x 1) containng values of the disturbances steps 
        %             used to gather step responses        
        %         
        % @note given step responses are onlly scaled by values and
        %       prolonged to the longest's response length. They are no 
        %       trimmed from the beggining and end. It is user's duty to
        %       do that
        %
        %================================================================
        function updateModel(obj, stepResponsesStruct)
            
            % Alias argument's name
            srs = stepResponsesStruct;
            
            % Set object's size
            obj.nu = size(srs.stepResponses, 1);
            obj.ny = size(srs.stepResponses, 2);
            
            % Establish value of the longest response
            D_t = 0;
            for u = 1:obj.nu
                for y = 1:obj.ny
                    if size(srs.stepResponses{u, y}, 1) > D_t
                        D_t = size(srs.stepResponses{u, y}, 1);
                    end
                end
            end
            
            % Establish scaling factors
            if ~isfield(srs, 'workPoint')
                srs.workPoint = zeros(obj.ny, 1);
            end
            if ~isfield(srs, 'stepSize')
                srs.stepSize = ones(obj.nu, 1);
            end
            
            % Extend step responses to a dynamic rank value
            obj.step_responses_shadow = zeros(obj.nu, obj.ny, D_t);
            for u = 1:obj.nu
                for y = 1:obj.ny
                    length = size(srs.stepResponses{u, y}, 1);
                    obj.step_responses_shadow(u, y, 1:length) = ...
                        (srs.stepResponses{u, y} - srs.workPoint(y)) / srs.stepSize(u);
                    obj.step_responses_shadow(u, y, length:end) = obj.step_responses_shadow(u, y, length);
                end
            end
            
            % If disturbance inputs present
            if isfield(srs, 'zStepResponses')
                
                % Set number of disturbance inputs
                obj.nz = size(srs.zStepResponses, 1);
            
                % Establish value of the longest response
                Dz_t = 0;
                for z = 1:obj.nz
                    for y = 1:obj.ny
                        if size(srs.zStepResponses{z, y}, 1) > Dz_t
                            Dz_t = size(srs.zStepResponses{z, y}, 1);
                        end
                    end
                end

                % Set step sizes
                if ~isfield(srs, 'zStepSize')
                    srs.zStepSize = ones(obj.nz, 1);
                end

                % Extend step responses to a dynamic rank value
                obj.z_step_responses_shadow = zeros(obj.nu, obj.ny, Dz_t);
                for z = 1:obj.nz
                    for y = 1:obj.ny
                        length = size(srs.zStepResponses{z, y}, 1);
                        obj.z_step_responses_shadow(z, y, 1:length) = ...
                            (srs.zStepResponses{z, y} - srs.workPoint(y)) / srs.zStepSize(z);
                        obj.z_step_responses_shadow(z, y, length:end) = obj.z_step_responses_shadow(z, y, length);
                    end
                end
                
            else
                
                obj.Dz_shadow = 0;
                obj.z_step_responses_shadow = [];
                
            end
            
            % Set regulator initialized 
            obj.initialized = true;
            
            % Reset internal
            obj.verify_params();
            obj.solve();
            obj.reset();
            
        end
        
        
        %================================================================
        % Function (re)initializes regulator using the function function
        % representing regulated object. It performs appropriate actions
        % to gether step responses and solve DMC optimisation problem
        % with them
        %
        % @param objectStruct : Structure containing information
        %        used to gather all required step responses. 
        %        The structure should containt following field:
        %
        %        .object - handle to the "function Y = f(U)"
        %                  or "function Y = f(U, Z) function representing 
        %                  regulated object. It should take two 2D arrays
        %                  representing CVs and DVs signals (matrices'
        %                  columns refer to a single CV/DV, when rows
        %                  contain signals valyes in subsequent discrete
        %                  moments). It should return output matrix of
        %                  size 'size(U, 1) x number of outputs' 
        %                  containing 'size(U, 1)' samples of outputs
        %
        %                  If number of rows of U and Z matrices are
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
        %                  values) that can be approximated as 0. This
        %                  value is used during step responses gathering
        %                  to establish whether the object's output is
        %                  stabilized
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
        % @note : If DMC parameters doesn't meet requirements described
        %         in the constructor's commet they are modified as
        %         described above (@see DMC(...))
        %
        % @note : Regulator's memory (i.e .uk_1, .dUp and .dZp vectors)
        %         are set to zero and should be updated after method's
        %         call if required. All limits are discarded.
        %================================================================
        function updateModelFunction(obj, objectStruct)
            
            % Compute step responses for CV-PV tracks
            step_responses_shadow_t   = cell(objectStruct.shape(1), objectStruct.shape(3));
            D_t = 0;
            for u = 1:objectStruct.shape(1)
                
                % Prepare SIMO bind
                bind_struct.object     = objectStruct.object;
                bind_struct.shape      = objectStruct.shape;
                bind_struct.U          = objectStruct.init_point_u;
                if isfield(objectStruct, 'init_point_z')
                    bind_struct.Z      = objectStruct.init_point_z;
                else
                    bind_struct.Z      = [];
                end                    
                bind_struct.input_num  = u;
                bind_struct.input_type = 'CV';

                % Initialize structure for single step response calculation
                srs_SIMO.object = @(input)(DMC.bind(bind_struct, input));
                srs_SIMO.init_point = objectStruct.init_point_u(u);
                srs_SIMO.step_size  = objectStruct.step_size_u(u);
                if isfield(objectStruct, 'tol')
                    srs_SIMO.tol = objectStruct.tol;
                end
                if isfield(objectStruct, 'sim_time')
                    srs_SIMO.sim_time = sobjectStructrs.sim_time;
                end

                % Compute step response
                step_responses_SIMO = DMC.step_response_SIMO(srs_SIMO);
                for y = 1:objectStruct.shape(3)
                    
                    step_responses_shadow_t{u, y} = step_responses_SIMO(:, y);
                    
                    % Update dynamic_rank
                    if size(step_responses_shadow_t{u, y}, 1) > D_t
                        D_t = size(step_responses_shadow_t{u, y}, 1);
                    end
                end

            end
            
            % Extend step responses to a dynamic rank value
            obj.step_responses_shadow = zeros(objectStruct.shape(1), objectStruct.shape(3), D_t);
            for u = 1:objectStruct.shape(1)
                for y = 1:objectStruct.shape(3)
                    length = size(step_responses_shadow_t{u, y}, 1);
                    obj.step_responses_shadow(u, y, 1:length) = step_responses_shadow_t{u, y};
                    obj.step_responses_shadow(u, y, length:end) = obj.step_responses_shadow(u, y, length);
                end
            end
            
            % Check if DVs are taken into account
            if objectStruct.shape(2) ~= 0

                % Compute step responses for DV-PV tracks
                z_step_responses_shadow_t = cell(objectStruct.shape(2), objectStruct.shape(3));
                Dz_t = 0;
                for z = 1:objectStruct.shape(2)

                    % Prepare SIMO bind
                    bind_struct.object     = objectStruct.object;
                    bind_struct.shape      = objectStruct.shape;
                    bind_struct.U          = objectStruct.init_point_u;
                    bind_struct.Z          = objectStruct.init_point_z;
                    bind_struct.input_num  = z;
                    bind_struct.output_num = y;
                    bind_struct.input_type = 'DV';

                    % Initialize structure for single step response calculation
                    srs_SIMO.object = @(input)(DMC.bind(bind_struct, input));
                    srs_SIMO.init_point = objectStruct.init_point_z(z);
                    srs_SIMO.step_size  = objectStruct.step_size_z(z);
                    if isfield(objectStruct, 'tol')
                        srs_SIMO.tol = objectStruct.tol;
                    end
                    if isfield(objectStruct, 'sim_time')
                        srs_SIMO.sim_time = objectStruct.sim_time;
                    end

                    % Compute step response
                    z_step_responses_SIMO = DMC.step_response_SIMO(srs_SIMO);
                    for y = 1:objectStruct.shape(3)

                        z_step_responses_shadow_t{z, y} = z_step_responses_SIMO(:, y);

                        % Update dynamic_rank
                        if size(z_step_responses_shadow_t{z, y}, 1) > Dz_t
                            Dz_t = size(z_step_responses_shadow_t{z, y}, 1);
                        end
                    end

                end
                
                % Extend step responses to a disturbance dynamic rank value
                obj.z_step_responses_shadow = zeros(objectStruct.shape(2), objectStruct.shape(3), Dz_t);
                for z = 1:objectStruct.shape(2)
                    for y = 1:objectStruct.shape(3)
                        length = size(z_step_responses_shadow_t{z, y}, 1);
                        obj.z_step_responses_shadow(z, y, 1:length) = z_step_responses_shadow_t{z, y};
                        obj.z_step_responses_shadow(z, y, length:end) = obj.z_step_responses_shadow(z, y, length);
                    end
                end

            else
                
                obj.Dz_shadow = 0;
                obj.z_step_responses_shadow = [];

            end
            
            % Set auxiliary values of object's shape
            obj.nu = size(obj.step_responses_shadow, 1);
            obj.nz = size(obj.z_step_responses_shadow, 1);
            obj.ny = size(obj.step_responses_shadow, 2);
            
            % Set regulator initialized 
            obj.initialized = true;
            
            % Reset internal
            obj.verify_params();
            obj.solve();
            obj.reset();
            
        end
        
        %================================================================
        % Computes optimal CV value for the given process' output and
        % desired output on the prediction horizon.
        %
        % @param y_zad : matrix of size (X x ny) representing desired 
        %        PVs values for next X discrete moments. Rows of the
        %        matrix refer to subsequent outputs when columns mark
        %        discrete moments.
        %        If X < N the Y_zad vector is extended to the
        %        X * ny length with the last element of the trajectory
        %        for each PV.
        %        If X > N the Y_zad vector is cut to X*ny.
        %
        % @param y : vector of shape (ny x 1) of actual PVs values
        %
        % @param z : vector of shape (nz x 1) of actual DVs values
        %================================================================
        function u = compute(obj, y_zad, y, z)
               
            % Check vector sized
            if size(y_zad, 2) ~= obj.ny
               error('y_zad vector is of the wrong size'); 
            end
            if size(y, 2) ~= obj.ny
               error('y vector is of the wrong size'); 
            end
            
            % If noise is taken into account update dZp property
            if obj.Dz ~= 0 
                if exist('z', 'var')

                    % Verify z vector size
                    if size(z, 1) ~= obj.nz
                       error('z vector is of the wrong size'); 
                    end

                    % Circshift dZp vector
                    obj.dZp = circshift(obj.dZp, obj.nz);

                    % Compute DVs  inreases
                    obj.dZp(1:obj.nz) = z - obj.zk_1;
                    
                end
            end
                
            % Initialize reshaped y_zad vector
            Y_zad = zeros(obj.N * obj.ny, 1);
            
            % Reshape and cut y_zad vector
            if size(y_zad, 1) > obj.N
                
                % Fill reshaped y_zad vector
                for i = 1:obj.N
                   Y_zad( ...
                       (i - 1) * obj.ny + 1 : ...
                            i  * obj.ny       ...
                    ) = y_zad(i, :)';
                       
                end

            % Reshape y_zad with or without extension
            else
                
                % Reshape y_zad
                for i = 1:size(y_zad, 1)
                   Y_zad( ...
                       (i - 1) * obj.ny + 1 : ...
                            i  * obj.ny       ...
                    ) = y_zad(i, :)';
                end
                
                % Extend y_zad
                if size(y_zad, 1) < obj.N
                    for i = size(y_zad, 1):obj.N
                        Y_zad( ...
                            (i - 1) * obj.ny + 1 : ...
                                 i  * obj.ny       ...
                         ) = y_zad(end, :)';
                    end
                end
            end
            
            % Validate y vector
            if size(y, 2) ~= obj.ny
               error('Too few output values given!') 
            end
            
            % Create Y vector
            Y = repmat(y', obj.N, 1);
        
            %---------------- Regulation law computation -----------------
            
            % Compute unrestrained trajectory
            if obj.Dz > 0
                Y_0 = Y + obj.Mp * obj.dUp + obj.Mzp * obj.dZp;
            else
                Y_0 = Y + obj.Mp * obj.dUp;
            end
            
            % Compute CVs increases
            du = obj.K(1:obj.nu, :) * (Y_zad - Y_0);
            
            %------------------------------------------------------------            
            
            % Apply dU value constraints
            for i = 1:size(du)
                if du(i) > obj.dumax(i)
                    du(i) = obj.dumax(i);
                elseif du(i) < obj.dumin(i)
                    du(i) = obj.dumin(i);
                end
            end

            % Circshift dUp vector
            obj.dUp = circshift(obj.dUp, obj.nu);

            % Save CVs  inreases
            obj.dUp(1:obj.nu) = du;
            
            % Compute CVs
            u = obj.uk_1 + du;
            
            % Apply U value constraints
            for i = 1:obj.nu
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
            obj.uk_1 = zeros(obj.nu, 1);
            obj.zk_1 = zeros(obj.nz, 1);
            obj.dUp  = zeros(obj.nu * (obj.D - 1), 1);
            obj.dZp  = zeros(obj.nz * obj.Dz     , 1);

            % Reset limits of controlled variables
            obj.umin   = ones(obj.nu, 1) * (-Inf);
            obj.umax   = ones(obj.nu, 1) * Inf;
            obj.dumin  = ones(obj.nu, 1) * (-Inf);
            obj.dumax  = ones(obj.nu, 1) * Inf;
            
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
            if obj.D > size(obj.step_responses_shadow, 3)
                obj.D_shadow =  size(obj.step_responses_shadow, 3);
                modified = true;
            end
            
            if obj.Dz > size(obj.z_step_responses_shadow, 3)
                obj.Dz_shadow = size(obj.z_step_responses_shadow, 3);
                modified = true;
            elseif obj.Dz > 0 && obj.nz == 0
                obj.Dz_shadow = 0;
                modified = true;
            end
            
            if obj.N > obj.D
                warning('N changed! Psi can be of the wrong size!')
                obj.N_shadow =  obj.D;
                modified = true;
            end
            
            if obj.Nu > obj.D
                warning('Nu changed! Lambda can be of the wrong size!')
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
                    obj.ny * obj.N, ...
                    obj.nu * obj.Nu ...
                 );
            for i = 1:obj.N
                for j = 1:obj.Nu
                    if i >= j
                        
                        % Create S matrix
                        S = obj.S(i - j + 1);
                        
                        % Fill field in M matrix
                        M((i-1)*obj.ny + 1 : i*obj.ny, ...
                          (j-1)*obj.nu + 1 : j*obj.nu) = S;
                      
                    end
                end
            end

            % Initialize Mp matrix
            obj.Mp = zeros( ...
                        obj.ny * obj.N,     ...
                        obj.nu *(obj.D - 1) ...
                     );
            for j = 1:(obj.D - 1)
                
                % Create Sj matrix
                Sj = obj.S(j);
                
                for i = 1:obj.N
                    
                    % Create S matrix
                    S = obj.S(i+j);
                    
                    % Fill Mp's field with Ss' difference
                    obj.Mp((i-1)*obj.ny + 1 : i*obj.ny, ...
                           (j-1)*obj.nu + 1 : j*obj.nu) = S - Sj;

                end
            end
            

            % Initialize Mzp matrix
            if obj.Dz ~= 0
                obj.Mzp = zeros(....
                              size(obj.step_responses_shadow,   2) * obj.N, ...
                              obj.nz * obj.Dz ...
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
                            obj.Mzp((i-1)*obj.ny + 1 : i*obj.ny, ...
                                    (j-1)*obj.nu + 1 : j*obj.nu) = S - Sj;
                        else
                            obj.Mzp((i-1)*obj.ny + 1 : i*obj.ny, ...
                                    (j-1)*obj.nu + 1 : j*obj.nu) = S;
                        end

                    end
               end
            else
                obj.Mzp = [];
            end
            
            % If 'lambda' and 'psi' parameters are given with scalars,
            % extend them to matrices
            if size(obj.lambda, 1) == 1
                lambda_t = eye(obj.Nu * obj.nu) * obj.lambda; 
            else
                lambda_t = obj.lambda;               
            end
            if size(obj.psi, 1) == 1
                psi_t = eye(obj.N * obj.ny) * obj.psi; 
            else
                psi_t = obj.psi;
            end
             
            % Solve optimization problem
            obj.K = (M' * psi_t * M + lambda_t)^(-1) * M' * psi_t;
        end
        
        %================================================================
        % Auxiliary function that constructs S_l matrix based on the
        % l index
        %
        % @param l : index of the samples
        %================================================================
        function S_l = S(obj, l)
            S_l = zeros(obj.ny, obj.nu);
            if l <= obj.D
                for y = 1:obj.ny
                    for u = 1:obj.nu
                        S_l(y, u) = obj.step_responses_shadow(u, y, l);
                    end
                end
            else
                for y = 1:obj.ny
                    for u = 1:obj.nu
                        S_l(y, u) = obj.step_responses_shadow(u, y, obj.D);
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
            Sz_l = zeros(size(obj.z_step_responses_shadow, 2), obj.nz);
            if l <= obj.Dz
                for y = 1:size(obj.z_step_responses_shadow, 2)
                    for z = 1:obj.nz
                        Sz_l(y, z) = obj.z_step_responses_shadow(z, y, l);
                    end
                end
            else
                for y = 1:size(obj.z_step_responses_shadow, 2)
                    for z = 1:obj.nz
                        Sz_l(y, z) = obj.z_step_responses_shadow(z, y, obj.Dz);
                    end
                end
            end
        end
        
    end
    
    methods (Access = private, Static = true)
        
        %================================================================
        % Computes step response suitable to use with DMC algorithm
        % for a single input - multiple output track.
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
        function Y = step_response_SIMO(step_response_struct)

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
            U = ones(sim_time, 1) * srs.init_point;
            Y = srs.object(U);

            % Get a stabilization time and work point
            stabilization_time = 0;
            for y = 1:size(Y, 2)
                stabilization_time_t = find(abs(Y(:, y) - Y(end, y)) < tol, 1, 'first');
                if isempty(stabilization_time_t)
                   error("Object's track could not be stabilized!") 
                end
                if stabilization_time_t > stabilization_time
                   stabilization_time = stabilization_time_t; 
                end
            end
            Ypp = Y(stabilization_time, :);

            % Measure step response
            U(stabilization_time:end) = step_response_struct.init_point + step_response_struct.step_size;
            Y = step_response_struct.object(U);

            % Scale step response
            for y = 1:size(Y, 2)
                Y(:, y) = (Y(:, y) - Ypp(y)) / step_response_struct.step_size;
            end
            
            % Trim step response from left and right
            if stabilization_time ~= size(Y, 1)
                Y = Y(stabilization_time+1 : end, :);
            else
                error('Object could not be stabilized!')
            end
            
            stabilization_time = 0;
            for y = 1:size(Y, 2)
                stabilization_time_t = find(abs(Y(:, y) - Y(end, y)) < tol, 1, 'first');
                if isempty(stabilization_time_t)
                   error("Object's track could not be stabilized!") 
                end
                if stabilization_time_t > stabilization_time
                   stabilization_time = stabilization_time_t; 
                end
            end
            
            % Return trimmed output
            Y = Y(1:stabilization_time, :);
            
        end

        %================================================================
        % Binds MIMO object function to a single SIMO track with given
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
        %        .input_type - type of the free input ('CV' or 'DV') 
        %
        % @param input : free input course
        %================================================================
        function output = bind(bind_struct, input)
            
            % Initialize input
            U = zeros(size(input, 1), bind_struct.shape(1));
            Z = zeros(size(input, 1), bind_struct.shape(2));
            
            % Fill input matrices
            for u = 1:size(U, 2)
                U(:, u) = bind_struct.U(u);
            end
            for z = 1:size(Z, 2)
                Z(:, z) = bind_struct.Z(z);
            end
            
            % Set free CV
            if strcmp(bind_struct.input_type,'CV')
                U(:, bind_struct.input_num) = input;
            elseif strcmp(bind_struct.input_type, 'DV')
                Z(:, bind_struct.input_num) = input;
            else
                error('Wrong input type!')
            end
            
            % Compute bind
            if bind_struct.shape(2) ~= 0
                output = bind_struct.object(U, Z);
            else
                output = bind_struct.object(U);
            end
            
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
            warning('N changed! Psi matrix can of be the wrong size')
            obj.N_shadow = N;
            obj.verify_params()
            obj.solve();
            obj.reset();
        end
        
        function val = get.Nu(obj)
            val = obj.Nu_shadow;
        end
        function set.Nu(obj, Nu)
            warning('Nu changed! Lambda matrix can of be the wrong size')
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
            
            % Check if lambda is of the right size
            if any(size(lambda) ~= [obj.Nu * obj.nu, obj.Nu * obj.nu]) && ...
               any(size(lambda) ~= [1 1])
               error('Lambda matrix should be either of (Nu * nu) x (Nu * nu) size or a scalar!') 
            end
            % If lambda is given with a matrix check if it is diagonal
            if any(size(lambda) == [obj.Nu * obj.nu, obj.Nu * obj.nu])
                if ~isdiag(lambda)
                   error('Lambda is not diagonal!') 
                end
            end
            obj.lambda_shadow = lambda;
            obj.verify_params();
            obj.solve();
            obj.reset();
        end

        function val = get.psi(obj)
            val = obj.psi_shadow;
        end
        function set.psi(obj, psi)
            % Check if psi is of the right size
            if any(size(psi) ~= [obj.N * obj.ny, obj.N * obj.ny]) && ...
               any(size(psi) ~= [1 1])
               error('Psi matrix should be either of (N * ny) x (N * ny) size or a scalar!') 
            end
            % If lambda is given with a matrix check if it is diagonal
            if any(size(psi) == [obj.N * obj.ny, obj.N * obj.ny])
                if ~isdiag(psi)
                   error('Psi is not diagonal!') 
                end
            end
            obj.psi_shadow = psi;
            obj.verify_params();
            obj.solve();
            obj.reset();
        end
        
        function K = get.K(obj)
            K = obj.K_shadow;
        end
        function set.K(obj, K)
            obj.K_shadow = K;
        end        
        
        function step_responses = get.step_responses(obj)
            step_responses = obj.step_responses_shadow;
        end
        function z_step_responses = get.z_step_responses (obj)
            z_step_responses  = obj.z_step_responses_shadow ;
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

