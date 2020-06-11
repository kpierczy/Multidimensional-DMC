%========================================================================
% File name   : FuzzyDMC.m
% Date        : 28th April 2020
% Version     : 0.1.0
% Author      : Krzysztof Pierczyk
% Description : Generic class implementing MIMO non-linear fuzzy DMC
%               regulator. Class takes into account both controller (CVs)
%               and uncontroller (DVs) inputs (internally uses DMC class,
%               @see DMC)
%
%========================================================================

classdef FuzzyDMC < handle

    properties (Dependent = true)
    
        % Limits of controlled variables
        umin  (:, 1) double {mustBeReal} = -Inf
        umax  (:, 1) double {mustBeReal} =  Inf
        dumin (:, 1) double {mustBeReal} = -Inf
        dumax (:, 1) double {mustBeReal} =  Inf
        
        % Last CVs values
        uk_1  (:, 1) double {mustBeReal}
        
    end
    
    properties (Access = public)
        
        % Sets-defining value
        sets_base  (1, :) char {mustBeMember( sets_base, { ...
                                    'input', ...
                                    'output'
                               })} = 'output'
        
        
        % Type of a membership function
        membership (1, :) char {mustBeMember( membership, { ...
                                    'gaussian', ...
                               })} = 'gaussian'
    end
    
    properties (Access = private)
        
        % Regulated object shape (number of CVs, DVs and PVS)
        shape (:, 1) uint32 {mustBeInteger, mustBeNonnegative}
        
        % Handle to the object-representing function
        object (1, 1)
        
        % Array of local regulators
        regulators cell     
        
    end
    
    properties (GetAccess = public, SetAccess = private)
       
        % Parameters of subsequent fuzzy sets
        c     (:, 1) double {mustBeReal} = []
        sigma (:, 1) double {mustBeReal} = []   
        
    end

    
    %====================================================================
    %%%%%%%%%%%%%%%%%%%%%%%%%% Public methods %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %====================================================================
    
    methods (Access = public)

        %================================================================
        % Fuzzy DMC regulator initialization. 
        % 
        % @param regulator_specific_struct : Columns array of structures 
        %        describing specific parameters of local regulators.
        %        This structure should contain following fields:
        %
        %        .sr_point - array conatining initial values of CVs
        %             used to gather step responses
        %
        %        .sr_step_size - array containing sizes of steps for
        %             all CVs used to gather step responses
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
        %
        %         There are also varian fields that are required only
        %         if objects meets some requirements. These are :
        %
        %        .sr_point_z (required if DVs are present) - array 
        %             containing initial values of DVs used to
        %             gather step responses
        %
        %        .sr_step_size_z (required if DVs are present) - array 
        %             containing sizes of steps for all DVs used 
        %             to gather step responses
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
        % @param common_struct : parameters common to all local regula-
        %        tors. This is structure that should contain following
        %        fields:
        %
        %         .shape - array conataining number of CVs, DVs and PVs
        %
        %         .object - handle to the function representing 
        %              regulated object. For requirements @see DMC
        %
        %         .membership_fun - name of the function used to
        %              compute values in fuzzy sets. Available functions
        %                  * 'gaussian'
        %
        %         .tol (optional, defaul = 0.0001) - maximum difference
        %               between two values (e.g. subsequent output
        %               values) that can be approximated as 0
        %
        %         .sim_time (optional, default = 1000) - length of the
        %               simulation used to get step responses
        %
        %===============================================================
        function obj = FuzzyDMC(regulator_specific_struct, common_struct)
            
            % Alias arguments names
            rss = regulator_specific_struct;
            cs  = common_struct;
            
            % Initialize structures for DMC cration
            sr_structures     = cell(size(rss, 1), 1); 
            params_structures = cell(size(rss, 1), 1);
            
            % Fill sr_structures
            for i = 1:size(sr_structures, 1)
                sr_structures{i, 1}.object = cs.object;
                sr_structures{i, 1}.shape = cs.shape;
                sr_structures{i, 1}.init_point_u = rss{i, 1}.sr_point;
                sr_structures{i, 1}.step_size_u = rss{i, 1}.sr_step_size;
                if cs.shape(2) ~= 0
                    sr_structures{i, 1}.init_point_z = rss{i, 1}.sr_point_z;
                    sr_structures{i, 1}.step_size_z = rss{i, 1}.sr_step_size_z;
                end
                if isfield(cs, 'tol')
                    sr_structures{i, 1}.tol = cs.tol;
                end
                if isfield(cs, 'sim_time')
                    sr_structures{i, 1}.sim_time = cs.sim_time;
                end                
            end
            
            % Fill params_structures
            for i = 1:size(params_structures, 1)
                params_structures{i, 1}.N = rss{i, 1}.N;
                params_structures{i, 1}.Nu = rss{i, 1}.Nu;
                params_structures{i, 1}.D = rss{i, 1}.D;
                params_structures{i, 1}.Dz = rss{i, 1}.Dz;
                params_structures{i, 1}.lambda = rss{i, 1}.lambda;
            end
            
            % Set type of membership function
            obj.membership = cs.membership_fun;
            
            % Save object shape
            obj.shape = cs.shape;
                        
            % Save parameters o fmembership function
            switch( obj.membership )
                case 'gaussian'
                    for i = 1:size(rss, 1)
                       obj.c(i) = rss{i, 1}.c; 
                    end
                    for i = 1:size(rss, 1)
                       obj.sigma(i) = rss{i, 1}.sigma; 
                    end
            end
            
            % Crate regulators
            obj.regulators = cell(size(rss, 1), 1); 
            for i = 1:size(rss, 1)
                obj.regulators{i, 1} = DMC(sr_structures{i, 1}, params_structures{i, 1});
            end

        end
        
        %================================================================
        % Computes optimal CV value for the given process' output and
        % desired output on the prediction horizon.
        %
        % @param y_zad : column vector of ((X * ny) x 1) shape 
        %        representing desired PVs values for next X discrete
        %        moments. 
        %        If X < N the Y_zad vector is extended to the
        %        X * ny length with the last element of the trajectory
        %        for each PV.
        %        If X > N the Y_zad vector is cut to X*ny.
        %
        % @param y : vector of actual PVs values
        %
        % @param z : vector of DVs (required if  number of DVs is > 0)
        %================================================================
        function u = compute(obj, y_zad, y, z)
            
            % Gathe local CVs
            cvs = zeros(obj.shape(1), size(obj.regulators, 1));
            for i = 1:size(cvs, 1)
                if exist('z', 'var')
                    cvs(:, i) = obj.regulators{i, 1}.compute(y_zad, y, z); 
                else
                    cvs(:, i) = obj.regulators{i, 1}.compute(y_zad, y); 
                end
            end
            
            % Compute membership weights
            memberships = zeros(obj.shape(1), size(obj.regulators, 1));
            for i = 1:size(memberships, 1)
                switch (obj.membership)
                    case 'gaussian'
                        if obj.sets_base == "input"
                            memberships(:, i) = gaussmf(obj.uk_1, [obj.sigma(i), obj.c(i)]);
                        else
                            memberships(:, i) = gaussmf(y       , [obj.sigma(i), obj.c(i)]);
                        end
                end                
            end
            
            % Sum membership weights
            wights_sum = sum(memberships, 2);
            
            % Compute CVs
            u = diag(cvs * memberships') ./ wights_sum;
            
            % Save actual CVs
            obj.uk_1 = u;            
             
        end

    end
    
    %====================================================================
    %%%%%%%%%%%%%%%%%%%%%%%% Getters & setters %%%%%%%%%%%%%%%%%%%%%%%%%%
    %====================================================================
    
    methods
               
        function val = get.umin(obj)
            val = obj.regulators{1}.umin;
        end
        function set.umin(obj, umin)
            for i=1:size(obj.regulators, 1)
                obj.regulators{i}.umin = umin;
            end
        end
        
        function val = get.umax(obj)
            val = obj.regulators{1}.umax;
        end
        function set.umax(obj, umax)
            for i=1:size(obj.regulators, 1)
                obj.regulators{i}.umax = umax;
            end
        end
        
        function val = get.dumin(obj)
            val = obj.regulators{1}.dumin;
        end
        function set.dumin(obj, dumin)
            for i=1:size(obj.regulators, 1)
                obj.regulators{i}.dumin = dumin;
            end
        end
        
        function val = get.dumax(obj)
            val = obj.regulators{1}.dumax;
        end
        function set.dumax(obj, dumax)
            for i=1:size(obj.regulators, 1)
                obj.regulators{i}.dumax = dumax;
            end
        end
        
        function val = get.uk_1(obj)
            val = obj.regulators{1}.uk_1;
        end
        function set.uk_1(obj, uk_1)
            for i=1:size(obj.regulators, 1)
                obj.regulators{i}.uk_1 = uk_1;
            end
        end        
    end
    
end

