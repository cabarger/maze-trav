//!
//! main.cpp
//!
//! CS 121.Bolden
//! Clang 16.0.6
//! Caleb Barger
//! 10/20/23 
//! Linux x86_64
//! barg8397@vandals.uidaho.edu
//!
//! Maze traversal program
//!

#include <chrono>
#include <stdio.h>
#include <stdlib.h>
#include <curses.h>

#include "cabarger_cs121_include.h"
#include "queue.h"

struct Coords {
  u8 row;
  u8 col;
};

QueuePrototype(Coords); 
QueueImpl(Coords);

/// Sample walk entry. Used within this programs pathfinding 
/// algorithm. See "mazeUpdateCellDistancesToTarget"
struct MazeSWE { 
  u16 distance; 
  Coords coords; 
};

QueuePrototype(MazeSWE);
QueueImpl(MazeSWE);

typedef u8 MazeCellChar;
enum {
  MazeCellChar_empty = '.',
  MazeCellChar_wall = '#',
  MazeCellChar_start = 'S',
  MazeCellChar_goal = 'G',
};

typedef u8 MazeCellType;
enum {
  MazeCellType_empty = 0,
  MazeCellType_wall,
  MazeCellType_start,
  MazeCellType_goal,
};

typedef i32 MazeCellColor; 
enum {
  MazeCellColor_queued = 1,
  MazeCellColor_visited,
  MazeCellColor_start,
  MazeCellColor_goal,
  MazeCellColor_cursor,
};

struct MazeCell {
  u16 queued : 1; 
  u16 visited : 1;
  u16 distance : 12;
  u16 type : 2;
};

struct Maze {
  u8 rows; 
  u8 cols;
  u16 start_index;

  MazeCell* cells; 
};

/// (row, col) deltas for N,S,E,W
const i8 cardinal_dir_deltas[][2] = {
  {-1,  0}, // N
  { 1,  0}, // S
  { 0,  1}, // E
  { 0, -1}  // W
};

inline MazeCellType mazeCellTypeFromChar(const MazeCellChar c) {
  switch(c) {
      case MazeCellChar_empty: return MazeCellType_empty;
      case MazeCellChar_wall: return MazeCellType_wall; 
      case MazeCellChar_start: return MazeCellType_start;
      case MazeCellChar_goal: return MazeCellType_goal;
      default: InvalidPath;
  }
  return 0; // NOTE(caleb): Shhhh -Wreturn-type 
}

inline MazeCellChar mazeCellCharFromType(const MazeCellType t) {
  switch(t) {
      case MazeCellType_empty: return MazeCellChar_empty;
      case MazeCellType_wall: return  MazeCellChar_wall; 
      case MazeCellType_start: return MazeCellChar_start;
      case MazeCellType_goal: return  MazeCellChar_goal;
      default: InvalidPath;
  }
  return 0; // NOTE(caleb): Shhhh -Wreturn-type 
}

inline bool validMazeCoords(Maze* m, i8 row, i8 col) {
  return ((row >= 0 && row < m->rows) && 
          (col >= 0 && col < m->cols)); 
}

inline MazeCell* mazeCellFromCoords(Maze* m, u8 row, u8 col) {
  const u16 cell_index = (u16)row * m->cols + col; 
  return (m->cells + cell_index);
}

Maze mazeFromStringU8(Arena* arena, StringU8 maze_str) {
  Maze result;
  u8 maze_validity = 0;

  u8 dims_written = 0;
  u16 maze_cell_index = 0;
  u16 current = 0;
  u16 start;
  while(current < maze_str.len) {
    start = current;
    const u8 c = maze_str.str[current++];
    switch(c) {
      case MazeCellChar_empty: 
      case MazeCellChar_wall: 
      case MazeCellChar_start:
      case MazeCellChar_goal: {
        const u8 row = maze_cell_index / result.cols; 
        const u8 col = maze_cell_index % result.cols;
        const MazeCellType cell_type = mazeCellTypeFromChar(c);
        result.cells[(u16)row * result.cols + col].type = cell_type;
        maze_validity |= (1 << cell_type);

        if (c == MazeCellChar_start) 
          result.start_index = maze_cell_index;

        maze_cell_index++;
      } break;
      
      case '0' ... '9': {
        while(isDigit(maze_str.str[current]))       
          current++;

        // NULL termd copy of current number for atoi lest undefined 
        // behavior occur.
        i32 dim;
        {
          ArenaState arena_state = arenaBegin(arena);

          StringU8 dim_z = StringU8{  
            .str = arenaAlloc(arena, u8, current - start + 1), 
            .len = 0,
            .cap = (u16)(current - start + 1),
          };
          for (u16 char_index=start; 
               char_index < start + current - start; 
               ++char_index)
            dim_z.str[dim_z.len++] = maze_str.str[char_index];
          dim_z.str[dim_z.len++] = 0;

          dim = atoi((char*)dim_z.str);

          arenaEnd(arena, arena_state); 
        }

        // Write row or col dim. If we have both dims allocate cells.
        if (dims_written == 0) {
          result.rows = (u8)dim;
        } else if (dims_written == 1) {
          result.cols = (u8)dim;
          result.cells = 
            arenaAlloc(arena, MazeCell, result.rows * result.cols);
          for (u16 cell_index=0; 
               cell_index < result.rows * result.cols; 
               ++cell_index)
            result.cells[cell_index] = MazeCell{ 
              .queued = 0, 
              .visited = 0, 
              .distance = 0, 
              .type = 0
            };
        } else 
          InvalidPath; 
        dims_written++;
      } break;
    }
  }

  // Assert valid maze
  if (!(maze_validity & (1 << MazeCellType_goal))) InvalidPath;
  if (!(maze_validity & (1 << MazeCellType_start))) InvalidPath;

  return result;
}

void mazeUpdateCellDistancesToTarget(
  Arena* scratch, 
  Maze* m, 
  Coords target
) {
  ArenaState scratch_restore = arenaBegin(scratch);

  // Reset cell distances
  for (u16 cell_index=0; cell_index < m->rows * m->cols; ++cell_index)
    m->cells[cell_index].distance = 0;

  QueueMazeSWE swe_queue = queueMazeSWEInit();
  queueMazeSWEEnqueue(&swe_queue, scratch, 
    MazeSWE{.distance = 0, .coords = target});
  while (swe_queue.head != 0) {
    MazeSWE swe = queueMazeSWEDequeue(&swe_queue);
    for (u8 dir_delta_index=0; 
         dir_delta_index < 4; 
         ++dir_delta_index) 
    {
      const i8 d_row = cardinal_dir_deltas[dir_delta_index][0];
      const i8 d_col = cardinal_dir_deltas[dir_delta_index][1];
      if (validMazeCoords(m, (i8)swe.coords.row + d_row, 
        (i8)swe.coords.col + d_col)) {
        Coords pos_to_queue = Coords{
          .row = (u8)(swe.coords.row + d_row), 
          .col = (u8)(swe.coords.col + d_col)
        };
        MazeCell* cell = mazeCellFromCoords(m, pos_to_queue.row, 
          pos_to_queue.col); 

        // Only queue if the distance still needs to be calculated.
        // Target has a distance of 0 so ignore it.
        if ((cell->distance == 0) && 
            ((pos_to_queue.row != target.row) || 
             (pos_to_queue.col != target.col))) 
        {
          // Set wall distance to the largest representable value. 
          // Otherwise update distance now and queue up the next sample
          // walk entry.
          if (cell->type == MazeCellType_wall) {   
            cell->distance = (1 << 12) - 1; 
          } else {
            cell->distance = swe.distance + 1;

            MazeSWE new_swe;  
            new_swe.distance = swe.distance + 1;
            new_swe.coords = pos_to_queue;
            queueMazeSWEEnqueue(&swe_queue, scratch, new_swe);
          }
        }
      }
    }
  }

  arenaEnd(scratch, scratch_restore);
}

u64 getTimeMS() {
  using namespace std::chrono;
  auto t = system_clock::now();
  auto since_epoch = t.time_since_epoch();
  return duration_cast<milliseconds>(since_epoch).count();
}

void drawCell(const MazeCell* cell, u8 row, u8 col) {
  if (cell->type == MazeCellType_goal)
    attron(COLOR_PAIR(MazeCellColor_goal));
  if (cell->queued != 0)
    attron(COLOR_PAIR(MazeCellColor_queued));
  if (cell->visited != 0)
    attron(COLOR_PAIR(MazeCellColor_visited));
  if (cell->type == MazeCellType_start)
    attron(COLOR_PAIR(MazeCellColor_start));
	mvaddch(row, col, mazeCellCharFromType(cell->type));
  attroff(A_COLOR);
}

void drawCursor(u8 row, u8 col) {
  attron(COLOR_PAIR(MazeCellColor_cursor));
	mvaddch(row, col, '@'); 
  attroff(A_COLOR);
}

int main(int argc, char** argv) {
  if (argc < 2) {
    // TODO(caleb): Usage notes
    return 0;
  }
  
  Arena scratch_arena = arenaInit(KB(2));
  Arena maze_arena = arenaInit(KB(4)); 
  
  // Curses init   
  initscr();
	cbreak();
	keypad(stdscr, TRUE);
  curs_set(0); 
  nodelay(stdscr, 1); // NOTE(caleb): Don't block on getch

  start_color();
  init_pair(MazeCellColor_queued, COLOR_MAGENTA, COLOR_BLACK);
  init_pair(MazeCellColor_visited, COLOR_CYAN, COLOR_BLACK);
  init_pair(MazeCellColor_start, COLOR_YELLOW, COLOR_BLACK);
  init_pair(MazeCellColor_goal, COLOR_RED, COLOR_BLACK);
  init_pair(MazeCellColor_cursor, COLOR_GREEN, COLOR_BLACK);

  // Load maze from file   
  Maze m; 
  {
    ArenaState scratch_state = arenaBegin(&scratch_arena); 
    FILE* mazef = fopen(argv[1], "r");

    // NOTE(caleb): The 40x40 maze is ~1.7kb 
    StringU8 maze_str = stringU8FromFile(&scratch_arena, mazef, KB(2));
    fclose(mazef);

    m = mazeFromStringU8(&maze_arena, maze_str); 
    arenaEnd(&scratch_arena, scratch_state);
  }
  
  QueueCoords coords_q = queueCoordsInit(); 

  // Place cursor/target pos(s) at maze entry point.
  Coords cursor_pos = {
    .row = (u8)(m.start_index / m.cols), 
    .col = (u8)(m.start_index % m.cols)
  };
  Coords target_pos = cursor_pos;

  u8 ticks_per_second = 10;
  u64 last_time = getTimeMS();

  bool done = false;
  while (!done) {
    if (getch() == KEY_UP) 
      ticks_per_second = min(60, ticks_per_second + 1);
    else if (getch() == KEY_DOWN) 
      ticks_per_second = max(1, ticks_per_second - 1);

    // Tick speed
    u64 now = getTimeMS();
    if (now - last_time < 1000 / ticks_per_second) continue;
    last_time = now;

    if ((cursor_pos.row == target_pos.row) && 
        (cursor_pos.col == target_pos.col)) 
    {
      // Since cursor has reached target, mark cell as visited
      MazeCell* curr_cell = mazeCellFromCoords(&m, 
        cursor_pos.row, cursor_pos.col);
      curr_cell->visited = 1; 
       
      // Add unvisited neighbors
      for (u8 dir_delta_index=0; dir_delta_index < 4; ++dir_delta_index) {
        const i8 d_row = cardinal_dir_deltas[dir_delta_index][0];
        const i8 d_col = cardinal_dir_deltas[dir_delta_index][1];
        if (validMazeCoords(&m, (i8)cursor_pos.row + d_row, 
          (i8)cursor_pos.col + d_col)) {
          Coords pos_to_queue = Coords{
            .row = (u8)(cursor_pos.row + d_row), 
            .col = (u8)(cursor_pos.col + d_col)
          };
          MazeCell* cell = mazeCellFromCoords(&m, 
            pos_to_queue.row, pos_to_queue.col); 
          if (cell->queued != 1) {
            if (cell->type != MazeCellType_wall)  {
              queueCoordsEnqueue(&coords_q, 
                &scratch_arena, pos_to_queue);
              cell->queued = 1;
            }
          }
        }
      }

      // Pull a new target off of the queue and update distances.
      target_pos = queueCoordsDequeue(&coords_q);
      mazeUpdateCellDistancesToTarget(&scratch_arena, &m, target_pos);
    } 

    // Either en'route to target or target_pos was just updated.
    if ((cursor_pos.row != target_pos.row) || 
        (cursor_pos.col != target_pos.col)) 
    {
      u8 nearest_cardinal_dir_delta_index;
      u16 nearest_distance = (1 << 12) - 1;
      for (u8 dir_delta_index=0; 
           dir_delta_index < 4; 
           ++dir_delta_index) 
      {
        const i8 d_row = cardinal_dir_deltas[dir_delta_index][0];
        const i8 d_col = cardinal_dir_deltas[dir_delta_index][1];
        if (validMazeCoords(&m, (i8)cursor_pos.row + d_row, 
          (i8)cursor_pos.col + d_col)) {
          Coords adj_cell_pos = Coords{
            .row = (u8)(cursor_pos.row + d_row), 
            .col = (u8)(cursor_pos.col + d_col)
          };
          MazeCell* cell = mazeCellFromCoords(&m, 
            adj_cell_pos.row, adj_cell_pos.col); 
          if (cell->distance < nearest_distance) {
            nearest_cardinal_dir_delta_index = dir_delta_index;
            nearest_distance = cell->distance;
          }
        }
      }
      cursor_pos.row += 
        cardinal_dir_deltas[nearest_cardinal_dir_delta_index][0];
      cursor_pos.col += 
        cardinal_dir_deltas[nearest_cardinal_dir_delta_index][1];

      // Check adjacent cells for goal
      for (u8 dir_delta_index=0; dir_delta_index < 4; ++dir_delta_index) {
        const i8 d_row = cardinal_dir_deltas[dir_delta_index][0];
        const i8 d_col = cardinal_dir_deltas[dir_delta_index][1];
        if (validMazeCoords(&m, (i8)cursor_pos.row + d_row, 
          (i8)cursor_pos.col + d_col)) {
          MazeCell* cell = mazeCellFromCoords(&m, 
            cursor_pos.row + d_row, cursor_pos.col + d_col);
          if (cell->type == MazeCellType_goal) 
            done = true; 
        }
      }
    }

    // Draw 
    for (u8 row_index=0; row_index < m.rows; ++row_index) {
      for (u8 col_index=0; col_index < m.cols; ++col_index) {
        MazeCell* cell = mazeCellFromCoords(&m, row_index, col_index);
        drawCell(cell, row_index, col_index);
      }
    }
    drawCursor(cursor_pos.row, cursor_pos.col);
  	refresh();
  }

  nodelay(stdscr, 0); // re-enable delay 
  getch(); // Wait for key before quiting 
  endwin();  // Deinit curses
}
