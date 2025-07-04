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
#include "dynamixel_sdk/group_bulk_write.h"


typedef struct
{
  uint8_t     id;
  uint16_t    data_end;
  uint16_t    start_address;
  uint16_t    data_length;
  uint8_t     *data;
}DataList;

typedef struct
{
  int         port_num;
  int         protocol_version;

  int         data_list_length;

  uint8_t     is_param_changed;

  uint16_t    param_length;

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
}

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

int groupBulkWrite(int port_num, int protocol_version)
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
  groupData[group_num].is_param_changed = False;
  groupData[group_num].param_length = 0;
  groupData[group_num].data_list = 0;

  groupBulkWriteClearParam(group_num);

  return group_num;
}

void groupBulkWriteMakeParam(int group_num)
{
  int data_num, idx, c;
  int port_num = groupData[group_num].port_num;

  if (groupData[group_num].protocol_version == 1)
    return;

  if (size(group_num) == 0)
    return;

  groupData[group_num].param_length = 0;

  idx = 0;
  for (data_num = 0; data_num < groupData[group_num].data_list_length; data_num++)
  {
    if (groupData[group_num].data_list[data_num].id == NOT_USED_ID)
      continue;

    groupData[group_num].param_length += 1 + 2 + 2 + groupData[group_num].data_list[data_num].data_length;

    packetData[port_num].data_write = (uint8_t*)realloc(packetData[port_num].data_write, groupData[group_num].param_length * sizeof(uint8_t));

    packetData[port_num].data_write[idx++] = groupData[group_num].data_list[data_num].id;
    packetData[port_num].data_write[idx++] = DXL_LOBYTE(groupData[group_num].data_list[data_num].start_address);
    packetData[port_num].data_write[idx++] = DXL_HIBYTE(groupData[group_num].data_list[data_num].start_address);
    packetData[port_num].data_write[idx++] = DXL_LOBYTE(groupData[group_num].data_list[data_num].data_length);
    packetData[port_num].data_write[idx++] = DXL_HIBYTE(groupData[group_num].data_list[data_num].data_length);

    for (c = 0; c < groupData[group_num].data_list[data_num].data_length; c++)
    {
      packetData[port_num].data_write[idx++] = groupData[group_num].data_list[data_num].data[c];
    }
  }
}

uint8_t groupBulkWriteAddParam(int group_num, uint8_t id, uint16_t start_address, uint16_t data_length, uint32_t data, uint16_t input_length)
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
    groupData[group_num].data_list[data_num].data_length = data_length;
    groupData[group_num].data_list[data_num].start_address = start_address;
    groupData[group_num].data_list[data_num].data = (uint8_t *)calloc(groupData[group_num].data_list[data_num].data_length, sizeof(uint8_t));
    groupData[group_num].data_list[data_num].data_end = 0;
  }
  else
  {
    if (groupData[group_num].data_list[data_num].data_end + input_length > groupData[group_num].data_list[data_num].data_length)
      return False;
  }

  switch (input_length)
  {
    case 1:
      groupData[group_num].data_list[data_num].data[groupData[group_num].data_list[data_num].data_end + 0] = DXL_LOBYTE(DXL_LOWORD(data));
      break;

    case 2:
      groupData[group_num].data_list[data_num].data[groupData[group_num].data_list[data_num].data_end + 0] = DXL_LOBYTE(DXL_LOWORD(data));
      groupData[group_num].data_list[data_num].data[groupData[group_num].data_list[data_num].data_end + 1] = DXL_HIBYTE(DXL_LOWORD(data));
      break;

    case 4:
      groupData[group_num].data_list[data_num].data[groupData[group_num].data_list[data_num].data_end + 0] = DXL_LOBYTE(DXL_LOWORD(data));
      groupData[group_num].data_list[data_num].data[groupData[group_num].data_list[data_num].data_end + 1] = DXL_HIBYTE(DXL_LOWORD(data));
      groupData[group_num].data_list[data_num].data[groupData[group_num].data_list[data_num].data_end + 2] = DXL_LOBYTE(DXL_HIWORD(data));
      groupData[group_num].data_list[data_num].data[groupData[group_num].data_list[data_num].data_end + 3] = DXL_HIBYTE(DXL_HIWORD(data));
      break;

    default:
      return False;
  }
  groupData[group_num].data_list[data_num].data_end = input_length;

  groupData[group_num].is_param_changed = True;
  return True;
}
void groupBulkWriteRemoveParam(int group_num, uint8_t id)
{
  int data_num = find(group_num, id);

  if (groupData[group_num].protocol_version == 1)
    return;

  if (data_num == groupData[group_num].data_list_length)
    return;

  if (groupData[group_num].data_list[data_num].id == NOT_USED_ID)  // NOT exist
    return;

  groupData[group_num].data_list[data_num].data_end = 0;

  groupData[group_num].data_list[data_num].data = 0;

  groupData[group_num].data_list[data_num].data_length = 0;
  groupData[group_num].data_list[data_num].start_address = 0;
  groupData[group_num].data_list[data_num].id = NOT_USED_ID;

  groupData[group_num].is_param_changed = True;
}

uint8_t groupBulkWriteChangeParam(int group_num, uint8_t id, uint16_t start_address, uint16_t data_length, uint32_t data, uint16_t input_length, uint16_t data_pos)
{
  int data_num = find(group_num, id);

  if (groupData[group_num].protocol_version == 1)
    return False;

  if (id == NOT_USED_ID)
    return False;

  if (data_num == groupData[group_num].data_list_length)
    return False;

  if (data_pos + input_length > groupData[group_num].data_list[data_num].data_length)
    return False;

  groupData[group_num].data_list[data_num].data_length = data_length;
  groupData[group_num].data_list[data_num].start_address = start_address;

  switch (input_length)
  {
    case 1:
      groupData[group_num].data_list[data_num].data[data_pos + 0] = DXL_LOBYTE(DXL_LOWORD(data));
      break;

    case 2:
      groupData[group_num].data_list[data_num].data[data_pos + 0] = DXL_LOBYTE(DXL_LOWORD(data));
      groupData[group_num].data_list[data_num].data[data_pos + 1] = DXL_HIBYTE(DXL_LOWORD(data));
      break;

    case 4:
      groupData[group_num].data_list[data_num].data[data_pos + 0] = DXL_LOBYTE(DXL_LOWORD(data));
      groupData[group_num].data_list[data_num].data[data_pos + 1] = DXL_HIBYTE(DXL_LOWORD(data));
      groupData[group_num].data_list[data_num].data[data_pos + 2] = DXL_LOBYTE(DXL_HIWORD(data));
      groupData[group_num].data_list[data_num].data[data_pos + 3] = DXL_HIBYTE(DXL_HIWORD(data));
      break;

    default:
      return False;
  }

  groupData[group_num].is_param_changed = True;
  return True;
}
void groupBulkWriteClearParam(int group_num)
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
void groupBulkWriteTxPacket(int group_num)
{
  if (groupData[group_num].protocol_version == 1)
  {
    packetData[groupData[group_num].port_num].communication_result = COMM_NOT_AVAILABLE;
    return;
  }

  if (size(group_num) == 0)
  {
    packetData[groupData[group_num].port_num].communication_result = COMM_NOT_AVAILABLE;
    return;
  }

  if (groupData[group_num].is_param_changed == True)
    groupBulkWriteMakeParam(group_num);

  bulkWriteTxOnly(groupData[group_num].port_num, groupData[group_num].protocol_version, groupData[group_num].param_length);
}
