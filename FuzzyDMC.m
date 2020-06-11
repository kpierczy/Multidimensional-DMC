%========================================================================
% File name   : FuzzyDMC.m
% Date        : 11th June 2020
% Version     : 0.1.0
% Author      : Krzysztof Pierczyk
% Description : Generic class implementing MIMO non-linear fuzzy DMC
%               regulator. Class takes into account both controller (CVs)
%               and uncontroller (DVs) inputs (internally uses DMC class,
%               @see DMC)
%
% @todo : implement other shapes of fuzzy sets
%
%========================================================================

classdef FuzzyDMC < handle

    properties (Access = public, Dependent = true)
        
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
    
    properties (Access = public)
        
        %================================================================
        % Sets-defining value. Decision about whether and how much 
        % actual state of the regulated object belongs to the fuzzy set
        % can be made basing on the input or the output from the object.
        %================================================================
        sets_base  (1, :) char {mustBeMember( sets_base, { ...
                                    'input', ...
                                    'output'
                               })} = 'output'
        
        
                           
        % Type of a membership function of fuuzzy sets
        membership (1, :) char {mustBeMember( membership, { ...
                                    'gaussian', ...
                               })} = 'gaussian'
                           
       % Array of local regulators
        regulators cell 
        
        % note : local regulators objects are public to enable easy
        %        parameters chaning. It should be emhpasized that
        %        local regulator's models should not be changed
        %        by the user after FuzzyDMC object's creation
        
    end
    
    properties (SetAccess = private, GetAccess = public)
        
        % Object's shape
        nu uint32 {mustBeNonnegative} = 0
        nz uint32 {mustBeNonnegative} = 0
        ny uint32 {mustBeNonnegative} = 0
       
        % Parameters of subsequent fuzzy sets
        c     (:, 1) double {mustBeReal} = []
        sigma (:, 1) double {mustBeReal} = []   
        
    end
    
    properties (Access = private)
        
        % Internals registers for imits of controlled variables
        umin  (:, 1) double {mustBeReal} = -Inf
        umax  (:, 1) double {mustBeReal} =  Inf
        dumin (:, 1) double {mustBeReal} = -Inf
        dumax (:, 1) double {mustBeReal} =  Inf
        
        % Last CVs values
        uk_1_shadow  (:, 1) double {mustBeReal}
        
        % Last DVs values
        zk_1_shadow  (:, 1) double {mustBeReal}
        
        % Vector of differences of subsequent CVs on the steerage horizon.
        dUp_shadow   (:, 1) double {mustBeReal} 
        
        % Vector of differences of subsequent DVs on the steerage horizon.
        dZp_shadow  (:, 1) double {mustBeReal} 
    
        D_top	(1, 1) uint32 {mustBeInteger, mustBePositive}    = 1
        Dz_top	(1, 1) uint32 {mustBeInteger, mustBeNonnegative} = 0
        
    end
        
    %====================================================================
    %%%%%%%%%%%%%%%%%%%%%%%%%% Public methods %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %====================================================================
    
    methods (Access = public)

        %================================================================
        % Fuzzy DMC regulator initialization. 
        % 
        % @param regulatorsParams : 1D cell array of structures 
        %        describing local DMCs' parameters. These structures
        %        should contain fields required by the DMC class'
        %        constructor (@see DMC)
        %
        % @param modelsStruct: 1D cell array of structures 
        %        describing initializing parameters of local regulators.
        %        Each structure should contain all fields required by
        %        either DMC.updateModel() or DMC.updateModelFunction
        %        (@see DMC) - an appropriate initializer is called
        %        depending on content of the particular structure.
        %        Additionally, each structure should contain some
        %        bonus fields:
        %
        %        .c (required if membership function is set to 
        %             'gaussian') - center of the fuzzy set that a local
        %             regulator belongs to
        %
        %        .sigma (required if membership function is set to 
        %             'gaussian') - standard deviation of the fuzzy set
        %             that a local regulator belongs to
        %
        %
        % @param membershipFun - name of the function used to
        %        compute values in fuzzy sets. Available functions:
        %           - 'gaussian'
        %
        %===============================================================
        function obj = FuzzyDMC(regulatorsParams, modelsStruct, membershipFun)
            
            % Check if number of regulator is established explicitly
            if size(regulatorsParams, 1) ~= size(modelsStruct)
               error('Sizes of regulatorsParams and modelsStruct are not equal!') 
            end
            
            % Alias arguments names
            ms = modelsStruct;
            
            % Set type of membership function
            if exist('membershipFun', 'var')
                obj.membership = membershipFun;
            end

            % Crate regulators
            obj.regulators = cell(size(ms, 1), 1); 
            for i = 1:size(ms, 1)
                obj.regulators{i} = DMC(regulatorsParams{i});
                if isfield(modelsStruct{i}, 'stepResponsesStruct')
                    obj.regulators{i}.updateModel(modelsStruct{i})
                else
                    obj.regulators{i}.updateModelFunction(modelsStruct{i})
                end
            end
            
            % Copy shape of the object
            obj.nu = obj.regulators{i}.nu;
            obj.nz = obj.regulators{i}.nz;
            obj.ny = obj.regulators{i}.ny;
            
            % Save parameters o fmembership function
            switch( obj.membership )
                case 'gaussian'
                    for i = 1:size(ms, 1)
                       obj.c(i) = ms{i, 1}.c; 
                    end
                    for i = 1:size(ms, 1)
                       obj.sigma(i) = ms{i, 1}.sigma; 
                    end
            end    
            
            % Find maximum values of D and Dz in the set of regulators
            obj.D_top = 1;
            for i = 1:size(obj.regulators, 1)
               if size(obj.regulators{i}.step_responses, 3) > obj.D_top
                   obj.D_top = size(obj.regulators{i}.step_responses, 3);
               end
            end
            obj.Dz_top = 0;
            if size(obj.regulators{1}.z_step_responses, 1) ~= 0
                for i = 1:size(obj.regulators, 1)
                   if size(obj.regulators{i}.z_step_responses, 3) > obj.Dz_top
                       obj.Dz_top = size(obj.regulators{i}.z_step_responses, 3);
                   end
                end
            end
            
            % Reset regulator's state
            obj.reset()

        end
        
        %================================================================
        % Computes optimal CV value for the given process' output and
        % desired output on the prediction horizon.
        %
        % @param y_zad : matrix of (X x ny) shape representing desired 
        %       PVs values for next X discrete moments. 
        %        If X < N the Y_zad vector is extended to the
        %        X * ny length with the last element of the trajectory
        %        for each PV.
        %        If X > N the Y_zad vector is cut to X*ny.
        %
        % @param y : vector of shape (ny x 1) of actual PVs values
        %
        % @param z : vector of shape (ny x 1) of actual DVs values
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
            if obj.Dz_top ~= 0 
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
            
            % Gather local CVs
            cvs = zeros(obj.nu, size(obj.regulators, 1));
            for i = 1:size(cvs, 1)
                if exist('z', 'var')
                    cvs(:, i) = obj.regulators{i, 1}.compute(y_zad, y, z); 
                else
                    cvs(:, i) = obj.regulators{i, 1}.compute(y_zad, y); 
                end
            end
            
            % Compute membership weights
            memberships = zeros(obj.nu, size(obj.regulators, 1));
            for i = 1:size(memberships, 1)
                switch (obj.membership)
                    case 'gaussian'
                        if obj.sets_base == "input"
                            memberships(:, i) = gaussmf(obj.uk_1, [obj.sigma(i), obj.c(i)]);
                        else
                            memberships(:, i) = gaussmf(       y, [obj.sigma(i), obj.c(i)]);
                        end
                end                
            end
            
            % Sum membership weights
            wights_sum = sum(memberships, 2);
            
            % Compute CVs
            u = sum(cvs .* memberships, 2) ./ wights_sum;
            
            % Apply U value constraints
            for i = 1:obj.nu
                if u(i) > obj.umax(i)
                    u(i) = obj.umax(i);
                elseif u(i) < obj.umin(i)
                    u(i) = obj.umin(i);
                end
            end
            
            % Compute CVs delta
            du = u - obj.uk_1;
            
            % Circshift dUp vector
            obj.dUp = circshift(obj.dUp, obj.nu);

            % Save CVs  inreases
            obj.dUp(1:obj.nu) = du;
            for i = 1:size(obj.regulators, 1)
                obj.regulators{i}.dUp(1:obj.nu) = du;
            end
            
            % Save actual CVs
            obj.uk_1 = u;    
            for i = 1:size(obj.regulators, 1)
                obj.regulators{i}.uk_1 = obj.uk_1;
            end
            if exist('z', 'var')
                obj.zk_1 = z;
                for i = 1:size(obj.regulators, 1)
                    obj.regulators{i}.zk_1 = obj.zk_1;
                end
            end 
            
        end

        %================================================================
        % Resets regulator's limits and memory to the default state.
        %================================================================
        function reset(obj)
           
            % Clear memory
            obj.uk_1_shadow = zeros(obj.nu, 1);
            obj.zk_1_shadow = zeros(obj.nz, 1);
            obj.dUp_shadow  = zeros(obj.nu * (obj.D_top - 1), 1);
            obj.dZp_shadow  = zeros(obj.nz * obj.Dz_top     , 1);

            % Reset limits of controlled variables
            obj.umin   = ones(obj.nu, 1) * (-Inf);
            obj.umax   = ones(obj.nu, 1) * Inf;
            obj.dumin  = ones(obj.nu, 1) * (-Inf);
            obj.dumax  = ones(obj.nu, 1) * Inf;
            
            % Reset subregulators
            for i = 1:size(obj.regulators, 1)
                obj.regulators{i}.reset()
            end
            
        end
        
    end
        
    %====================================================================
    %%%%%%%%%%%%%%%%%%%%%%%% Getters & setters %%%%%%%%%%%%%%%%%%%%%%%%%%
    %====================================================================
    
    methods
        
        function val = get.uk_1(obj)
            val = obj.uk_1_shadow;
        end
        function set.uk_1(obj, uk_1)
            obj.uk_1_shadow = uk_1;
            for i=1:size(obj.regulators, 1)
                obj.regulators{i}.uk_1 = obj.uk_1_shadow;
            end
        end   
        
        function val = get.zk_1(obj)
            val = obj.zk_1_shadow;
        end
        function set.zk_1(obj, zk_1)
            obj.zk_1_shadow = zk_1;
            for i=1:size(obj.regulators, 1)
                obj.regulators{i}.zk_1 = obj.zk_1_shadow;
            end
        end   
        
        function val = get.dUp(obj)
            val = obj.dUp_shadow;
        end
        function set.dUp(obj, dUp)
            obj.dUp_shadow = dUp;
            for i=1:size(obj.regulators, 1)
                obj.regulators{i}.dUp = ...
                    obj.dUp_shadow(1:obj.nu * (obj.regulators{i}.D - 1));
            end
        end   
        
        function val = get.dZp(obj)
            val = obj.dUp_shadow;
        end
        function set.dZp(obj, dZp)
            obj.dZp_shadow = dZp;
            for i=1:size(obj.regulators, 1)
                obj.regulators{i}.dZp = ...
                    obj.dZp_shadow(obj.nz * obj.regulators{i}.Dz);
            end
        end
    end
    
end

