/*******************************************************************************
* Copyright (c) 2016, ROBOTIS CO., LTD.
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
*
* * Redistributions of source code must retain the above copyright notice, this
*   list of conditions and the following disclaimer.
*
* * Redistributions in binary form must reproduce the above copyright notice,
*   this list of conditions and the following disclaimer in the documentation
*   and/or other materials provided with the distribution.
*
* * Neither the name of ROBOTIS nor the names of its
*   contributors may be used to endorse or promote products derived from
*   this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
* FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
* SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
* CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
* OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
* OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*******************************************************************************/

/* Author: Ryu Woon Jung (Leon) */

#if defined(_WIN32) || defined(_WIN64)
#define WINDLLEXPORT
#endif

#include <stdlib.h>
#include "dynamixel_sdk/group_sync_read.h"


typedef struct
{
  uint8_t     id;
  uint8_t     *data;
}DataList;

typedef struct
{
  int         port_num;
  int         protocol_version;

  int         data_list_length;

  uint8_t     last_result;
  uint8_t     is_param_changed;

  uint16_t    start_address;
  uint16_t    data_length;

  DataList   *data_list;
}GroupData;

static GroupData *groupData;
static int g_used_group_num = 0;

static int size(int group_num)
{
  int data_num;
  int real_size = 0;

  for (data_num = 0; data_num < groupData[group_num].data_list_length; data_num++)
  {
    if (groupData[group_num].data_list[data_num].id != NOT_USED_ID)
      real_size++;
  }
  return real_size;
};

static int find(int group_num, int id)
{
  int data_num;

  for (data_num = 0; data_num < groupData[group_num].data_list_length; data_num++)
  {
    if (groupData[group_num].data_list[data_num].id == id)
      break;
  }

  return data_num;
}

int groupSyncRead(int port_num, int protocol_version, uint16_t start_address, uint16_t data_length)
{
  int group_num = 0;

  if (g_used_group_num != 0)
  {
    for (group_num = 0; group_num < g_used_group_num; group_num++)
    {
      if (groupData[group_num].is_param_changed != True)
        break;
    }
  }

  if (group_num == g_used_group_num)
  {
    g_used_group_num++;
    groupData = (GroupData *)realloc(groupData, g_used_group_num * sizeof(GroupData));
  }

  groupData[group_num].port_num = port_num;
  groupData[group_num].protocol_version = protocol_version;
  groupData[group_num].data_list_length = 0;
  groupData[group_num].last_result = False;
  groupData[group_num].is_param_changed = False;
  groupData[group_num].start_address = start_address;
  groupData[group_num].data_length = data_length;
  groupData[group_num].data_list = 0;

  groupSyncReadClearParam(group_num);

  return group_num;
}

void groupSyncReadMakeParam(int group_num)
{
  int data_num, idx;
  int port_num = groupData[group_num].port_num;

  if (groupData[group_num].protocol_version == 1)
    return;

  if (size(group_num) == 0)
    return;

  packetData[port_num].data_write = (uint8_t*)realloc(packetData[port_num].data_write, size(group_num) * (1) * sizeof(uint8_t)); // ID(1)

  idx = 0;
  for (data_num = 0; data_num < groupData[group_num].data_list_length; data_num++)
  {
    if (groupData[group_num].data_list[data_num].id == NOT_USED_ID)
      continue;

    packetData[port_num].data_write[idx++] = groupData[group_num].data_list[data_num].id;
  }
}

uint8_t groupSyncReadAddParam(int group_num, uint8_t id)
{
  int data_num = 0;

  if (groupData[group_num].protocol_version == 1)
    return False;

  if (id == NOT_USED_ID)
    return False;

  if (groupData[group_num].data_list_length != 0)
    data_num = find(group_num, id);

  if (groupData[group_num].data_list_length == data_num)
  {
    groupData[group_num].data_list_length++;
    groupData[group_num].data_list = (DataList *)realloc(groupData[group_num].data_list, groupData[group_num].data_list_length * sizeof(DataList));

    groupData[group_num].data_list[data_num].id = id;
    groupData[group_num].data_list[data_num].data = (uint8_t *)calloc(groupData[group_num].data_length, sizeof(uint8_t));
  }

  groupData[group_num].is_param_changed = True;
  return True;
}
void groupSyncReadRemoveParam(int group_num, uint8_t id)
{
  int data_num = find(group_num, id);

  if (groupData[group_num].protocol_version == 1)
    return;

  if (groupData[group_num].data_list[data_num].id == NOT_USED_ID)  // NOT exist
    return;

  groupData[group_num].data_list[data_num].data = 0;

  groupData[group_num].data_list[data_num].id = NOT_USED_ID;

  groupData[group_num].is_param_changed = True;
}
void groupSyncReadClearParam(int group_num)
{
  int port_num = groupData[group_num].port_num;

  if (groupData[group_num].protocol_version == 1)
    return;

  if (size(group_num) == 0)
    return;

  groupData[group_num].data_list = 0;

  packetData[port_num].data_write = 0;

  groupData[group_num].data_list_length = 0;

  groupData[group_num].is_param_changed = False;
}

void groupSyncReadTxPacket(int group_num)
{
  int port_num = groupData[group_num].port_num;

  if (groupData[group_num].protocol_version == 1)
  {
    packetData[port_num].communication_result = COMM_NOT_AVAILABLE;
    return;
  }

  if (size(group_num) == 0)
  {
    packetData[port_num].communication_result = COMM_NOT_AVAILABLE;
    return;
  }

  if (groupData[group_num].is_param_changed == True)
    groupSyncReadMakeParam(group_num);

  syncReadTx(groupData[group_num].port_num
    , groupData[group_num].protocol_version
    , groupData[group_num].start_address
    , groupData[group_num].data_length
    , (size(group_num) * 1));
}

void groupSyncReadRxPacket(int group_num)
{
  int data_num, c;
  int port_num = groupData[group_num].port_num;

  groupData[group_num].last_result = False;

  if (groupData[group_num].protocol_version == 1)
  {
    packetData[port_num].communication_result = COMM_NOT_AVAILABLE;
    return;
  }

  packetData[groupData[group_num].port_num].communication_result = COMM_RX_FAIL;

  if (size(group_num) == 0)
  {
    packetData[groupData[group_num].port_num].communication_result = COMM_NOT_AVAILABLE;
    return;
  }

  for (data_num = 0; data_num < groupData[group_num].data_list_length; data_num++)
  {
    if (groupData[group_num].data_list[data_num].id == NOT_USED_ID)
      continue;

      packetData[port_num].data_read
        = (uint8_t *)realloc(packetData[port_num].data_read, groupData[group_num].data_length * sizeof(uint8_t));

      readRx(groupData[group_num].port_num, groupData[group_num].protocol_version, groupData[group_num].data_length);
      if (packetData[port_num].communication_result != COMM_SUCCESS)
        return;

      for (c = 0; c < groupData[group_num].data_length; c++)
        groupData[group_num].data_list[data_num].data[c] = packetData[port_num].data_read[c];
  }

  if (packetData[port_num].communication_result == COMM_SUCCESS)
    groupData[group_num].last_result = True;
}

void groupSyncReadTxRxPacket(int group_num)
{
  int port_num = groupData[group_num].port_num;

  if (groupData[group_num].protocol_version == 1)
  {
    packetData[port_num].communication_result = COMM_NOT_AVAILABLE;
    return;
  }

  packetData[port_num].communication_result = COMM_TX_FAIL;

  groupSyncReadTxPacket(group_num);
  if (packetData[port_num].communication_result != COMM_SUCCESS)
    return;

  groupSyncReadRxPacket(group_num);
}

uint8_t groupSyncReadIsAvailable(int group_num, uint8_t id, uint16_t address, uint16_t data_length)
{
  int data_num = find(group_num, id);

  if (groupData[group_num].protocol_version == 1 || groupData[group_num].last_result == False || groupData[group_num].data_list[data_num].id == NOT_USED_ID)
    return False;

  if (address < groupData[group_num].start_address || groupData[group_num].start_address + groupData[group_num].data_length - data_length < address) {
    return False;
  }
  return True;
}

uint32_t groupSyncReadGetData(int group_num, uint8_t id, uint16_t address, uint16_t data_length)
{
  int data_num = find(group_num, id);

  if (groupSyncReadIsAvailable(group_num, id, address, data_length) == False)
    return 0;

  switch (data_length)
  {
    case 1:
      return groupData[group_num].data_list[data_num].data[address - groupData[group_num].start_address];

    case 2:
      return DXL_MAKEWORD(groupData[group_num].data_list[data_num].data[address - groupData[group_num].start_address], groupData[group_num].data_list[data_num].data[address - groupData[group_num].start_address + 1]);

    case 4:
      return DXL_MAKEDWORD(
        DXL_MAKEWORD(
        groupData[group_num].data_list[data_num].data[address - groupData[group_num].start_address + 0]
        , groupData[group_num].data_list[data_num].data[address - groupData[group_num].start_address + 1])
        , DXL_MAKEWORD(
          groupData[group_num].data_list[data_num].data[address - groupData[group_num].start_address + 2]
          , groupData[group_num].data_list[data_num].data[address - groupData[group_num].start_address + 3])
      );

    default:
      return 0;
  }
}
