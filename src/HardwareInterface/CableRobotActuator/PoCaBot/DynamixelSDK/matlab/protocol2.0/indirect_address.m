%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2016, ROBOTIS CO., LTD.
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% * Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.
%
% * Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution.
%
% * Neither the name of ROBOTIS nor the names of its
%   contributors may be used to endorse or promote products derived from
%   this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Author: Ryu Woon Jung (Leon)

%
% *********     Indirect Address Example      *********
%
%
% Available Dynamixel model on this example : All models using Protocol 2.0
% This example is designed for using a Dynamixel PRO 54-200, and an USB2DYNAMIXEL.
% To use another Dynamixel model, such as X series, see their details in E-Manual(support.robotis.com) and edit below variables yourself.
% Be sure that Dynamixel PRO properties are already set as %% ID : 1 / Baudnum : 3 (Baudrate : 1000000 [1M])
%

clc;
clear all;

lib_name = '';

if strcmp(computer, 'PCWIN')
  lib_name = 'dxl_x86_c';
elseif strcmp(computer, 'PCWIN64')
  lib_name = 'dxl_x64_c';
elseif strcmp(computer, 'GLNX86')
  lib_name = 'libdxl_x86_c';
elseif strcmp(computer, 'GLNXA64')
  lib_name = 'libdxl_x64_c';
end

% Load Libraries
if ~libisloaded(lib_name)
    [notfound, warnings] = loadlibrary(lib_name, 'dynamixel_sdk.h', 'addheader', 'port_handler.h', 'addheader', 'packet_handler.h', 'addheader', 'group_sync_write.h', 'addheader', 'group_sync_read.h');
end

% Control table address                                            % Control table address is different in Dynamixel model
ADDR_PRO_INDIRECTADDRESS_FOR_WRITE      = 49;                  % EEPROM region
ADDR_PRO_INDIRECTADDRESS_FOR_READ       = 59;                  % EEPROM region
ADDR_PRO_TORQUE_ENABLE                  = 562;
ADDR_PRO_LED_RED                        = 563;
ADDR_PRO_GOAL_POSITION                  = 596;
ADDR_PRO_MOVING                         = 610;
ADDR_PRO_PRESENT_POSITION               = 611;
ADDR_PRO_INDIRECTDATA_FOR_WRITE         = 634;
ADDR_PRO_INDIRECTDATA_FOR_READ          = 639;

% Data Byte Length
LEN_PRO_LED_RED                         = 1;
LEN_PRO_GOAL_POSITION                   = 4;
LEN_PRO_MOVING                          = 1;
LEN_PRO_PRESENT_POSITION                = 4;
LEN_PRO_INDIRECTDATA_FOR_WRITE          = 5;
LEN_PRO_INDIRECTDATA_FOR_READ           = 5;

% Protocol version
PROTOCOL_VERSION                        = 2.0;                % See which protocol version is used in the Dynamixel

% Default setting
DXL_ID                                  = 1;                  % Dynamixel ID: 1
BAUDRATE                                = 1000000;
DEVICENAME                              = 'COM1';             % Check which port is being used on your controller
                                                              % ex) Windows: "COM1"   Linux: "/dev/ttyUSB0"

TORQUE_ENABLE                           = 1;                  % Value for enabling the torque
TORQUE_DISABLE                          = 0;                  % Value for disabling the torque
DXL_MINIMUM_POSITION_VALUE              = -150000;            % Dynamixel will rotate between this value
DXL_MAXIMUM_POSITION_VALUE              = 150000;             % and this value (note that the Dynamixel would not move when the position value is out of movable range. Check e-manual about the range of the Dynamixel you use.)
DXL_MINIMUM_LED_VALUE                   = 0;                  % Dynamixel LED will light between this value
DXL_MAXIMUM_LED_VALUE                   = 255;                % and this value
DXL_MOVING_STATUS_THRESHOLD             = 20;                 % Dynamixel moving status threshold

ESC_CHARACTER                           = 'e';                % Key for escaping loop

COMM_SUCCESS                            = 0;                  % Communication Success result value
COMM_TX_FAIL                            = -1001;              % Communication Tx Failed

% Initialize PortHandler Structs
% Set the port path
% Get methods and members of PortHandlerLinux or PortHandlerWindows
port_num = portHandler(DEVICENAME);

% Initialize PacketHandler Structs
packetHandler();

% Initialize Groupsyncwrite Structs
groupwrite_num = groupSyncWrite(port_num, PROTOCOL_VERSION, ADDR_PRO_INDIRECTDATA_FOR_WRITE, LEN_PRO_INDIRECTDATA_FOR_WRITE);

% Initialize Groupsyncread Structs for Present Position
groupread_num = groupSyncRead(port_num, PROTOCOL_VERSION, ADDR_PRO_INDIRECTDATA_FOR_READ, LEN_PRO_INDIRECTDATA_FOR_READ);

index = 1;
dxl_comm_result = COMM_TX_FAIL;                               % Communication result
dxl_addparam_result = false;                                  % AddParam result
dxl_getdata_result = false;                                   % GetParam result
dxl_goal_position = [DXL_MINIMUM_POSITION_VALUE DXL_MAXIMUM_POSITION_VALUE];         % Goal position

dxl_error = 0;                                                % Dynamixel error
dxl_moving = 0;                                               % Dynamixel moving status
dxl_led_value = [0 255];                                      % Dynamixel LED value
dxl_present_position = 0;                                     % Present position


% Open port
if (openPort(port_num))
    fprintf('Succeeded to open the port!\n');
else
    unloadlibrary(lib_name);
    fprintf('Failed to open the port!\n');
    input('Press any key to terminate...\n');
    return;
end


% Set port baudrate
if (setBaudRate(port_num, BAUDRATE))
    fprintf('Succeeded to change the baudrate!\n');
else
    unloadlibrary(lib_name);
    fprintf('Failed to change the baudrate!\n');
    input('Press any key to terminate...\n');
    return;
end


% Disable Dynamixel Torque :
% Indirect address would not accessible when the torque is already enabled
write1ByteTxRx(port_num, PROTOCOL_VERSION, DXL_ID, ADDR_PRO_TORQUE_ENABLE, TORQUE_DISABLE);
if getLastTxRxResult(port_num, PROTOCOL_VERSION) ~= COMM_SUCCESS
    printTxRxResult(PROTOCOL_VERSION, getLastTxRxResult(port_num, PROTOCOL_VERSION));
elseif getLastRxPacketError(port_num, PROTOCOL_VERSION) ~= 0
    printRxPacketError(PROTOCOL_VERSION, getLastRxPacketError(port_num, PROTOCOL_VERSION));
else
    fprintf('Dynamixel has been successfully connected \n');
end

% INDIRECTDATA parameter storages replace LED, goal position, present position and moving status storages
write2ByteTxRx(port_num, PROTOCOL_VERSION, DXL_ID, ADDR_PRO_INDIRECTADDRESS_FOR_WRITE + 0, ADDR_PRO_GOAL_POSITION + 0);
if getLastTxRxResult(port_num, PROTOCOL_VERSION) ~= COMM_SUCCESS
    printTxRxResult(PROTOCOL_VERSION, getLastTxRxResult(port_num, PROTOCOL_VERSION));
elseif getLastRxPacketError(port_num, PROTOCOL_VERSION) ~= 0
    printRxPacketError(PROTOCOL_VERSION, getLastRxPacketError(port_num, PROTOCOL_VERSION));
end

write2ByteTxRx(port_num, PROTOCOL_VERSION, DXL_ID, ADDR_PRO_INDIRECTADDRESS_FOR_WRITE + 2, ADDR_PRO_GOAL_POSITION + 1);
if getLastTxRxResult(port_num, PROTOCOL_VERSION) ~= COMM_SUCCESS
    printTxRxResult(PROTOCOL_VERSION, getLastTxRxResult(port_num, PROTOCOL_VERSION));
elseif getLastRxPacketError(port_num, PROTOCOL_VERSION) ~= 0
    printRxPacketError(PROTOCOL_VERSION, getLastRxPacketError(port_num, PROTOCOL_VERSION));
end

write2ByteTxRx(port_num, PROTOCOL_VERSION, DXL_ID, ADDR_PRO_INDIRECTADDRESS_FOR_WRITE + 4, ADDR_PRO_GOAL_POSITION + 2);
if getLastTxRxResult(port_num, PROTOCOL_VERSION) ~= COMM_SUCCESS
    printTxRxResult(PROTOCOL_VERSION, getLastTxRxResult(port_num, PROTOCOL_VERSION));
elseif getLastRxPacketError(port_num, PROTOCOL_VERSION) ~= 0
    printRxPacketError(PROTOCOL_VERSION, getLastRxPacketError(port_num, PROTOCOL_VERSION));
end

write2ByteTxRx(port_num, PROTOCOL_VERSION, DXL_ID, ADDR_PRO_INDIRECTADDRESS_FOR_WRITE + 6, ADDR_PRO_GOAL_POSITION + 3);
if getLastTxRxResult(port_num, PROTOCOL_VERSION) ~= COMM_SUCCESS
    printTxRxResult(PROTOCOL_VERSION, getLastTxRxResult(port_num, PROTOCOL_VERSION));
elseif getLastRxPacketError(port_num, PROTOCOL_VERSION) ~= 0
    printRxPacketError(PROTOCOL_VERSION, getLastRxPacketError(port_num, PROTOCOL_VERSION));
end

write2ByteTxRx(port_num, PROTOCOL_VERSION, DXL_ID, ADDR_PRO_INDIRECTADDRESS_FOR_WRITE + 8, ADDR_PRO_LED_RED);
if getLastTxRxResult(port_num, PROTOCOL_VERSION) ~= COMM_SUCCESS
    printTxRxResult(PROTOCOL_VERSION, getLastTxRxResult(port_num, PROTOCOL_VERSION));
elseif getLastRxPacketError(port_num, PROTOCOL_VERSION) ~= 0
    printRxPacketError(PROTOCOL_VERSION, getLastRxPacketError(port_num, PROTOCOL_VERSION));
end

write2ByteTxRx(port_num, PROTOCOL_VERSION, DXL_ID, ADDR_PRO_INDIRECTADDRESS_FOR_READ + 0, ADDR_PRO_PRESENT_POSITION + 0);
if getLastTxRxResult(port_num, PROTOCOL_VERSION) ~= COMM_SUCCESS
    printTxRxResult(PROTOCOL_VERSION, getLastTxRxResult(port_num, PROTOCOL_VERSION));
elseif getLastRxPacketError(port_num, PROTOCOL_VERSION) ~= 0
    printRxPacketError(PROTOCOL_VERSION, getLastRxPacketError(port_num, PROTOCOL_VERSION));
end

write2ByteTxRx(port_num, PROTOCOL_VERSION, DXL_ID, ADDR_PRO_INDIRECTADDRESS_FOR_READ + 2, ADDR_PRO_PRESENT_POSITION + 1);
if getLastTxRxResult(port_num, PROTOCOL_VERSION) ~= COMM_SUCCESS
    printTxRxResult(PROTOCOL_VERSION, getLastTxRxResult(port_num, PROTOCOL_VERSION));
elseif getLastRxPacketError(port_num, PROTOCOL_VERSION) ~= 0
    printRxPacketError(PROTOCOL_VERSION, getLastRxPacketError(port_num, PROTOCOL_VERSION));
end

write2ByteTxRx(port_num, PROTOCOL_VERSION, DXL_ID, ADDR_PRO_INDIRECTADDRESS_FOR_READ + 4, ADDR_PRO_PRESENT_POSITION + 2);
if getLastTxRxResult(port_num, PROTOCOL_VERSION) ~= COMM_SUCCESS
    printTxRxResult(PROTOCOL_VERSION, getLastTxRxResult(port_num, PROTOCOL_VERSION));
elseif getLastRxPacketError(port_num, PROTOCOL_VERSION) ~= 0
    printRxPacketError(PROTOCOL_VERSION, getLastRxPacketError(port_num, PROTOCOL_VERSION));
end

write2ByteTxRx(port_num, PROTOCOL_VERSION, DXL_ID, ADDR_PRO_INDIRECTADDRESS_FOR_READ + 6, ADDR_PRO_PRESENT_POSITION + 3);
if getLastTxRxResult(port_num, PROTOCOL_VERSION) ~= COMM_SUCCESS
    printTxRxResult(PROTOCOL_VERSION, getLastTxRxResult(port_num, PROTOCOL_VERSION));
elseif getLastRxPacketError(port_num, PROTOCOL_VERSION) ~= 0
    printRxPacketError(PROTOCOL_VERSION, getLastRxPacketError(port_num, PROTOCOL_VERSION));
end

write2ByteTxRx(port_num, PROTOCOL_VERSION, DXL_ID, ADDR_PRO_INDIRECTADDRESS_FOR_READ + 8, ADDR_PRO_MOVING);
if getLastTxRxResult(port_num, PROTOCOL_VERSION) ~= COMM_SUCCESS
    printTxRxResult(PROTOCOL_VERSION, getLastTxRxResult(port_num, PROTOCOL_VERSION));
elseif getLastRxPacketError(port_num, PROTOCOL_VERSION) ~= 0
    printRxPacketError(PROTOCOL_VERSION, getLastRxPacketError(port_num, PROTOCOL_VERSION));
end

% Enable Dynamixel Torque
write1ByteTxRx(port_num, PROTOCOL_VERSION, DXL_ID, ADDR_PRO_TORQUE_ENABLE, TORQUE_ENABLE);
if getLastTxRxResult(port_num, PROTOCOL_VERSION) ~= COMM_SUCCESS
    printTxRxResult(PROTOCOL_VERSION, getLastTxRxResult(port_num, PROTOCOL_VERSION));
elseif getLastRxPacketError(port_num, PROTOCOL_VERSION) ~= 0
    printRxPacketError(PROTOCOL_VERSION, getLastRxPacketError(port_num, PROTOCOL_VERSION));
end

% Add parameter storage for the present position value
dxl_addparam_result = groupSyncReadAddParam(groupread_num, DXL_ID);
if dxl_addparam_result ~= true
  fprintf('[ID:%03d] groupSyncRead addparam failed', DXL1_ID);
  return;
end


while 1
    if input('Press any key to continue! (or input e to quit!)\n', 's') == ESC_CHARACTER
        break;
    end

    % Add values to the Syncwrite storage
    dxl_addparam_result = groupSyncWriteAddParam(groupwrite_num, DXL_ID, typecast(int32(dxl_goal_position(index)), 'uint32'), LEN_PRO_GOAL_POSITION);
    if dxl_addparam_result ~= true
        fprintf('[ID:%03d] groupSyncWrite addparam failed', DXL1_ID);
        return;
    end
    dxl_addparam_result = groupSyncWriteAddParam(groupwrite_num, DXL_ID, dxl_led_value(index), LEN_PRO_LED_RED);
    if dxl_addparam_result ~= true
        fprintf('[ID:%03d] groupSyncWrite addparam failed', DXL1_ID);
        return;
    end

    % Syncwrite all
    groupSyncWriteTxPacket(groupwrite_num);
    if getLastTxRxResult(port_num, PROTOCOL_VERSION) ~= COMM_SUCCESS
        printTxRxResult(PROTOCOL_VERSION, getLastTxRxResult(port_num, PROTOCOL_VERSION));
    end

    % Clear syncwrite parameter storage
    groupSyncWriteClearParam(groupwrite_num);

    while 1
      % Syncread present position from indirectdata2
      groupSyncReadTxRxPacket(groupread_num);
      if getLastTxRxResult(port_num, PROTOCOL_VERSION) ~= COMM_SUCCESS
          printTxRxResult(PROTOCOL_VERSION, getLastTxRxResult(port_num, PROTOCOL_VERSION));
      end

      % Check if groupsyncread data of Dyanamixel is available
      dxl_getdata_result = groupSyncReadIsAvailable(groupread_num, DXL_ID, ADDR_PRO_INDIRECTDATA_FOR_READ, LEN_PRO_PRESENT_POSITION);
      if dxl_getdata_result ~= true
        fprintf('[ID:%03d] groupSyncRead getdata failed', DXL1_ID);
        return;
      end

      % Check if groupsyncread data of Dyanamixel is available
      dxl_getdata_result = groupSyncReadIsAvailable(groupread_num, DXL_ID, ADDR_PRO_INDIRECTDATA_FOR_READ + LEN_PRO_PRESENT_POSITION, LEN_PRO_MOVING);
      if dxl_getdata_result ~= true
        fprintf('[ID:%03d] groupSyncRead getdata failed', DXL1_ID);
        return;
      end

      % Get Dynamixel present position value
      dxl_present_position = groupSyncReadGetData(groupread_num, DXL_ID, ADDR_PRO_INDIRECTDATA_FOR_READ, LEN_PRO_PRESENT_POSITION);

      % Get Dynamixel moving status value
      dxl_moving = groupSyncReadGetData(groupread_num, DXL_ID, ADDR_PRO_INDIRECTDATA_FOR_READ + LEN_PRO_PRESENT_POSITION, LEN_PRO_MOVING);

      fprintf('[ID:%03d] GoalPos:%d  PresPos:%d  IsMoving:%d\n', DXL_ID, dxl_goal_position(index), typecast(uint32(dxl_present_position), 'int32'), dxl_moving);

      if ~(abs(dxl_goal_position(index) - typecast(uint32(dxl_present_position), 'int32')) > DXL_MOVING_STATUS_THRESHOLD)
        break;
      end
    end

    % Change goal position
    if index == 1
        index = 2;
    else
        index = 1;
    end
end


% Disable Dynamixel Torque
write1ByteTxRx(port_num, PROTOCOL_VERSION, DXL_ID, ADDR_PRO_TORQUE_ENABLE, TORQUE_DISABLE);
if getLastTxRxResult(port_num, PROTOCOL_VERSION) ~= COMM_SUCCESS
    printTxRxResult(PROTOCOL_VERSION, getLastTxRxResult(port_num, PROTOCOL_VERSION));
elseif getLastRxPacketError(port_num, PROTOCOL_VERSION) ~= 0
    printRxPacketError(PROTOCOL_VERSION, getLastRxPacketError(port_num, PROTOCOL_VERSION));
end

% Close port
closePort(port_num);

% Unload Library
unloadlibrary(lib_name);

close all;
clear all;
