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

#include "cabarger_cs121_include.h"
#include "queue.h"

#include <curses.h>
#include <stdlib.h> 

// Coords/MazeSWE queue prototype(s)/implementation(s). 
// Via scuffed templates.

struct Coords {
  u8 row;
  u8 col;
};

QueuePrototype(Coords); 
QueueImpl(Coords);

/// Sample walk entry. Relavent to pathfinding algorithm used.  
struct MazeSWE { 
  u16 distance; 
  Coords coords; 
};

QueuePrototype(MazeSWE);
QueueImpl(MazeSWE);

////////////////////

typedef u8 MazeCellType;
enum {
  MazeCellType_empty = 0,
  MazeCellType_wall,
  MazeCellType_start,
  MazeCellType_goal,
};

/// Accompanying char representations for the 4 cell types. 
/// Only 2 bits per cell to encode all cell types.
typedef u8 MazeCellChar;
enum {
  MazeCellChar_empty = '.',
  MazeCellChar_wall = '#',
  MazeCellChar_start = 'S',
  MazeCellChar_goal = 'G',
};

// MazeCellChar <=> MazeCellType conversions

inline MazeCellType cellTypeFromChar(const MazeCellChar c) {
  switch(c) {
      case MazeCellChar_empty: return MazeCellType_empty;
      case MazeCellChar_wall: return MazeCellType_wall; 
      case MazeCellChar_start: return MazeCellType_start;
      case MazeCellChar_goal: return MazeCellType_goal;
      default: InvalidPath;
  }
  return 0; // NOTE(caleb): Shhhh -Wreturn-type 
}

inline MazeCellChar cellCharFromType(const MazeCellType t) {
  switch(t) {
      case MazeCellType_empty: return MazeCellChar_empty;
      case MazeCellType_wall: return  MazeCellChar_wall; 
      case MazeCellType_start: return MazeCellChar_start;
      case MazeCellType_goal: return  MazeCellChar_goal;
      default: InvalidPath;
  }
  return 0; // NOTE(caleb): Shhhh -Wreturn-type 
}

/// 2 byte bitfield encoding a cell. 
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
  MazeCell* data; 
};

struct AdjCellCoords {
  Coords coords[4];
  u8 count;  
};

/// Various forground/background pair ids for cell coloring.
typedef i32 MazeCellColor; 
enum {
  MazeCellColor_queued = 1,
  MazeCellColor_visited,
  MazeCellColor_start,
  MazeCellColor_goal,
  MazeCellColor_cursor,
};

////////////////////

// Maze operations

inline bool mazeBoundsCheck(Maze* m, i8 row, i8 col) {
  return ((row >= 0 && row < m->rows) && 
          (col >= 0 && col < m->cols)); 
}

inline MazeCell* mazeCellFromCoords(Maze* m, u8 row, u8 col) {
  const u16 cell_index = (u16)row * m->cols + col; 
  return (m->data + cell_index);
}

inline MazeCell* mazeCellFromCoords(Maze* m, Coords coords) {
  return mazeCellFromCoords(m, coords.row, coords.col);
}

AdjCellCoords mazeAdjCellCoords(Maze* m, Coords pos) {
  AdjCellCoords result;
  result.count = 0;

  // (row, col) deltas for N,S,E,W
  const i8 dir_deltas[][2] = {
    {-1,  0}, // N
    { 1,  0}, // S
    { 0,  1}, // E
    { 0, -1}  // W
  };

  for (u8 dir_delta_index=0; dir_delta_index < 4; ++dir_delta_index) {
    const i8 d_row = dir_deltas[dir_delta_index][0];
    const i8 d_col = dir_deltas[dir_delta_index][1];
    if (mazeBoundsCheck(m, (i8)pos.row + d_row, 
        (i8)pos.col + d_col)) {
      result.coords[result.count++] = Coords{
        .row = (u8)(pos.row + d_row), 
        .col = (u8)(pos.col + d_col)
      };
    }
  }
  
  return result;
}

Maze mazeFromStringU8(Arena* arena, StringU8 maze_str) {
  Maze result;

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
        result.data[(u16)row * result.cols + col].type =
           cellTypeFromChar(c);

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

          arenaEnd(arena, arena_state); // End dim_z scratch
        }

        // Write row or col and if we have both dims allocate 
        // space for cells.
        if (dims_written == 0) {
          result.rows = (u8)dim;
        } else if (dims_written == 1) {
          result.cols = (u8)dim;
          result.data = 
            arenaAlloc(arena, MazeCell, result.rows * result.cols);
        } else 
          InvalidPath; 
        dims_written++;
      } break;
    }
  }
  
  return result;
}

void mazeUpdateCellDistancesToTarget(
  Maze* m, 
  Arena* scratch, 
  Coords target
) {
  ArenaState scratch_restore = arenaBegin(scratch);

  for (u16 cell_index=0; cell_index < m->rows * m->cols; ++cell_index)
    m->data[cell_index].distance = 0;

  QueueMazeSWE swe_queue = queueMazeSWEInit();
  queueMazeSWEEnqueue(&swe_queue, scratch, 
    MazeSWE{.distance = 0, .coords = target});
  while (swe_queue.head != 0) {
    MazeSWE swe = queueMazeSWEDequeue(&swe_queue);

    AdjCellCoords adj_cell_coords = mazeAdjCellCoords(m, swe.coords);
    for (u8 adj_index=0; adj_index < adj_cell_coords.count; 
         ++adj_index) {
      MazeCell* cell = mazeCellFromCoords(m, 
        adj_cell_coords.coords[adj_index]); 

      // Only queue if the distance still needs to be calculated.
      // Target has a distance of 0 so ignore it.
      if ((cell->distance == 0) && 
          ((adj_cell_coords.coords[adj_index].row != target.row) || 
           (adj_cell_coords.coords[adj_index].col != target.col))) 
      {
        // Set wall distance to the largest representable value. 
        // Otherwise update distance now queue up the next sample
        // walk entry.
        if (cell->type == MazeCellType_wall) {   
          cell->distance = (1 << 12) - 1; 
        } else {
          cell->distance = swe.distance + 1;
          queueMazeSWEEnqueue(&swe_queue, scratch, MazeSWE{
            .distance = cell->distance,
            .coords = adj_cell_coords.coords[adj_index]
          });
        }
      }
    }
  }

  arenaEnd(scratch, scratch_restore);
}

////////////////////

// Drawing primitives

void drawCell(const MazeCell* cell, u8 row, u8 col) {
  if (cell->type == MazeCellType_goal)
    attron(COLOR_PAIR(MazeCellColor_goal));
  if (cell->queued != 0)
    attron(COLOR_PAIR(MazeCellColor_queued));
  if (cell->visited != 0)
    attron(COLOR_PAIR(MazeCellColor_visited));
  if (cell->type == MazeCellType_start)
    attron(COLOR_PAIR(MazeCellColor_start));
	mvaddch(row, col, cellCharFromType(cell->type));
  attroff(A_COLOR);
}

void drawCursor(u8 row, u8 col) {
  attron(COLOR_PAIR(MazeCellColor_cursor));
	mvaddch(row, col, '@'); 
  attroff(A_COLOR);
}

////////////////////

int main(int argc, char** argv) {
  if (argc < 2) {
    printf("Usage: %s path_to_maze\n", argv[0]);
    return 0;
  }
  
  Arena scratch_arena = arenaInit(KB(2));
  Arena maze_arena = arenaInit(KB(4)); 
  
  // Curses init   
  initscr();
	cbreak();
	keypad(stdscr, TRUE);
  curs_set(0); 
  nodelay(stdscr, 1);

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

  f32 ticks_per_second = 10.0;
  f64 last_time = getTimeMS();

  bool done = false;
  while (!done) {
    if (getch() == KEY_UP) 
      ticks_per_second = Min(60.0, ticks_per_second + 1.0);
    else if (getch() == KEY_DOWN) 
      ticks_per_second = Max(1.0, ticks_per_second - 1.0);

    // Tick speed 
    f64 now = getTimeMS();
    if (now - last_time < 1000.0 / ticks_per_second) continue;
    last_time = now;

    // Queue adj cells after reaching target pos.
    if ((cursor_pos.row == target_pos.row) && 
        (cursor_pos.col == target_pos.col)) 
    {
      // Since cursor has reached target, mark cell as visited. 
      // NOTE(caleb): "visited" field is ONLY used for coloring.
      {
        MazeCell* curr_cell = mazeCellFromCoords(&m, cursor_pos);
        curr_cell->visited = 1; 
      }
       
      // Queue adj cells that haven't yet been queued (and aren't walls)
      AdjCellCoords adj_cell_coords = mazeAdjCellCoords(&m, cursor_pos);
      for (u8 adj_index=0; adj_index < adj_cell_coords.count; 
           ++adj_index) {
        MazeCell* cell = mazeCellFromCoords(&m, 
          adj_cell_coords.coords[adj_index]);
        if (cell->queued != 1 && cell->type != MazeCellType_wall) {
          queueCoordsEnqueue(&coords_q, 
            &scratch_arena, adj_cell_coords.coords[adj_index]);
          cell->queued = 1;
        }
      }

      // Pull a new target off queue and update distances.
      target_pos = queueCoordsDequeue(&coords_q);
      mazeUpdateCellDistancesToTarget(&m, &scratch_arena, target_pos);
    } 

    // Either en'route to target or target_pos was just updated.
    if ((cursor_pos.row != target_pos.row) || 
        (cursor_pos.col != target_pos.col)) 
    {
      u8 nearest_adj_index = 0;
      u16 nearest_distance = (1 << 12) - 1;

      AdjCellCoords adj_cell_coords = mazeAdjCellCoords(&m, cursor_pos);
      for (u8 adj_index=0; adj_index < adj_cell_coords.count; 
           ++adj_index) {
        MazeCell* cell = mazeCellFromCoords(&m, 
          adj_cell_coords.coords[adj_index]);
        if (cell->distance < nearest_distance) {
          nearest_adj_index = adj_index;
          nearest_distance = cell->distance;
        }
      }
      cursor_pos = adj_cell_coords.coords[nearest_adj_index];

      // Check adjacent cells for goal
      adj_cell_coords = mazeAdjCellCoords(&m, cursor_pos);
      for (u8 adj_index=0; adj_index < adj_cell_coords.count; 
           ++adj_index) {
        MazeCell* cell = mazeCellFromCoords(&m, 
          adj_cell_coords.coords[adj_index]);
        if (cell->type == MazeCellType_goal) 
          done = true; 
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
